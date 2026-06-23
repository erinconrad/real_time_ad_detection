function [state, newADs, newStim] = rt_ad_detector_step(state, block)
% RT_AD_DETECTOR_STEP  Feed one block of raw EEG to the streaming detector.
%
%   [state, newADs, newStim] = rt_ad_detector_step(state, block)
%
% INPUT
%   state : from rt_ad_detector_init (and threaded through prior steps)
%   block : [nSamples x nChannels] new raw samples (microvolts), where
%           nChannels matches numel(chLabels). nSamples is arbitrary, so
%           you can deliver data in 100-1000 ms packets to emulate a true
%           real-time feed.
%
% OUTPUT
%   state   : updated state
%   newADs  : table (Type/Channels/OnTime/OffTime) of AD detections produced
%             by this block (OnTime = absolute s; OffTime = NaN, online)
%   newStim : table of stim onset rows produced by this block
%
% The detector advances internally in chunkDuration steps; the residual tail
% that doesn't fill a full analysis window is carried to the next block, so
% block boundaries do not affect the result (results are identical to
% feeding the whole session at once).

p   = state.p;
fs  = state.fs;
numCh = state.numChannels;

% append incoming samples to residual
res = [state.residual; block];
nRes = size(res,1);

newAD_idx   = [];
newStim_idx = [];

ptr = 1;
while ptr + state.chunkSize - 1 <= nRes
    window = res(ptr:ptr+state.chunkSize-1, :);
    global_end_index = state.samp_offset + ptr + state.chunkSize - 1; % 1-based
    abs_time = state.t0 + global_end_index/fs;

    [state, addedAD, addedStim] = process_window(state, window, abs_time);
    newAD_idx   = [newAD_idx;   addedAD];   %#ok<AGROW>
    newStim_idx = [newStim_idx; addedStim]; %#ok<AGROW>

    ptr = ptr + state.updateSize;
end

% carry remaining tail
consumed = ptr - 1;
state.residual    = res(consumed+1:end, :);
state.samp_offset = state.samp_offset + consumed;

newADs  = state.events(newAD_idx, :);
newStim = state.events(newStim_idx, :);

end


function [state, addedAD, addedStim] = process_window(state, dataChunk, abs_time)
% One internal analysis window. Mirrors the inner loop of find_ad_fcn.
p     = state.p;
fs    = state.fs;
numCh = state.numChannels;
addedAD   = [];
addedStim = [];

% periodic decay of bad-channel counter
state.n_bad_loop_counter = state.n_bad_loop_counter + 1;
if state.n_bad_loop_counter == p.n_bad_reset
    state.bad_ch_counter = max([zeros(1,numCh); ...
        state.bad_ch_counter - p.n_reduce_bad*ones(1,numCh)], [], 1);
    state.n_bad_loop_counter = 0;
end

% nan handling + demean
if all(isnan(dataChunk),'all'), dataChunk = zeros(size(dataChunk)); end
for ich = 1:numCh
    dataChunk(isnan(dataChunk(:,ich)),ich) = mean(dataChunk(:,ich),'omitnan');
end
dataChunk = dataChunk - mean(dataChunk,1,'omitnan');

% causal high-pass (stateful)
[dataChunk, state.zi] = stevefilter(dataChunk, state.zi, p.hpf_alpha);

% drop excluded channels (columnwise)
dataChunk(:,state.exclude) = nan;

% bipolar montage
alt   = state.altBipolarIndices;
valid = ~isnan(alt);
newDataChunk = dataChunk;
newDataChunk(:,valid)  = dataChunk(:,valid) - dataChunk(:,alt(valid));
newDataChunk(:,~valid) = nan;
newDataChunk = newDataChunk - mean(newDataChunk,1,'omitnan');

% decide look flags
look_for_off = 0;
if state.stim_on == 0
    look_for_stim = 1;
    if isempty(state.last_stim_chs)
        look_for_ad = 0;
    elseif abs_time > state.last_stim_off + p.stopLookingADSecs
        look_for_ad = 0; state.look_for_sat = 0;
    else
        look_for_ad = 1;
    end
else
    look_for_stim = 0; look_for_off = 1; look_for_ad = 0; state.look_for_sat = 1;
end

% decaying buffer (emphasizes repetitive 50 Hz stim)
state.buffer(isnan(state.buffer)) = dataChunk(isnan(state.buffer));
state.buffer = state.buffer*p.decay + dataChunk;
buffer_power = sum(state.buffer.^2,1);

%% Stim onset detection
if look_for_stim
    state.bad_ch_counter = state.bad_ch_counter + sum(newDataChunk > p.bad_ch_amp,1);

    chs_above = buffer_power > p.stimPowerBoost;
    state.last_ones_stim(1:end-1,:) = state.last_ones_stim(2:end,:);
    state.last_ones_stim(end,:) = chs_above;
    detected_stim = sum(state.last_ones_stim==1,1) > ...
        size(state.last_ones_stim,1)*p.perc_above_thresh_stim;

    [~,bip] = find_bipolar_pairs(state.chLabels(detected_stim), find(detected_stim));
    if ~isempty(bip)
        if size(bip,1) > 1
            mp = arrayfun(@(j) mean(buffer_power(bip(j,:))), 1:size(bip,1));
            [~,h] = max(mp);
            keep = bip(h,:);
        else
            keep = bip;
        end

        r1 = add_event(state, 'stim', state.chLabels{keep(1)}, abs_time);
        state = r1.state; row1 = r1.row;
        r2 = add_event(state, 'stim', state.chLabels{keep(2)}, abs_time);
        state = r2.state; row2 = r2.row;
        addedStim = [addedStim; row1; row2];

        state.last_ones_stim(:) = 0;
        state.stim_on        = 1;
        state.keep_pair      = keep;
        state.last_stim_on   = abs_time;
        state.last_stim_chs  = keep;
        state.last_stim_rows = [row1 row2];
        state.ad_chs_this_stim = [];
        state.baselines = state.baselines_all(1:p.n_baseline_keep);
    end

    % roll baseline history (use pre-stim windows)
    state.baselines_all(1:end-1) = state.baselines_all(2:end);
    state.baselines_all(end)     = {newDataChunk};
end

%% Stim offset detection
if look_for_off
    if mean(buffer_power(state.keep_pair)) < p.stimOffPower
        state.stim_on = 0;
        state.events.OffTime(state.last_stim_rows) = abs_time;
        state.last_stim_off = abs_time;
    end
end

%% Power features for AD detection
if state.look_for_sat || look_for_ad
    power = measure_power(newDataChunk);

    bl = cellfun(@(x) measure_power(x), state.baselines, 'UniformOutput', false);
    bl = mean(cell2mat(bl), 1);
    rel = power ./ bl;

    above = rel > p.ad_thresh & rel < p.ad_too_high_thresh;

    hf = measure_power(newDataChunk, p.hfband, fs);
    above(hf > p.hfthresh) = 0;

    saturated = rel > p.ad_too_high_thresh;
    state.time_since_last_sat(saturated) = abs_time;
    above(abs_time - state.time_since_last_sat < p.coolDownLastSat) = 0;
end

%% AD detection
if look_for_ad
    state.last_ones(1:end-1,:) = state.last_ones(2:end,:);
    state.last_ones(end,:) = above;

    detected_ad = sum(state.last_ones==1,1) > p.num_above_thresh;

    chs_to_look = ad_look_chs(state.chLabels, state.altBipolarIndices, ...
        state.last_stim_chs, state.file_name);
    detected_ad(chs_to_look == 0) = 0;
    detected_ad(state.bad_ch_counter > p.n_bad) = 0;

    if any(detected_ad)
        ones_idx = find(detected_ad);
        for k = 1:numel(ones_idx)
            if ismember(ones_idx(k), state.ad_chs_this_stim), continue; end
            r = add_event(state, 'AD', state.chLabels{ones_idx(k)}, abs_time);
            state = r.state;
            addedAD = [addedAD; r.row]; %#ok<AGROW>
            state.ad_chs_this_stim(end+1) = ones_idx(k);
        end
        if ~isempty(state.ad_chs_this_stim)
            state.last_ones(:, state.ad_chs_this_stim) = 0;
        end
    end
else
    state.last_ones = zeros(state.num_thresh, numCh);
end

end


function out = add_event(state, type, ch, onTime)
% Append a row to the master event log; return new state + row index.
state.events = [state.events; {type, ch, onTime, NaN}];
out.state = state;
out.row   = height(state.events);
end
