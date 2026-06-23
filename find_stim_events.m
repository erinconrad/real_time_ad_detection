function [stimT, params] = find_stim_events(session, params)
% FIND_STIM_EVENTS  Detect 50 Hz stimulation events (channel pair + on/off
% times) in a loaded session. Each detected event becomes one "trial" for
% the ground-truth annotation and the sensitivity/specificity analysis.
%
%   stimT = find_stim_events(session)
%   stimT = find_stim_events(session, params)
%
% INPUT
%   session : struct from export_sessions (.values, .fs, .chLabels,
%             .start_time, .ieeg_name)
%   params  : optional struct overriding defaults (see below)
%
% OUTPUT
%   stimT   : table with one row per stim event:
%       EventNum, StimChs {1x2 labels}, StimIdx [1x2], OnTime, OffTime (abs s)
%
% The detection logic mirrors the stim portion of the original find_ad_fcn:
% a decaying buffer emphasizes repetitive (50 Hz) signals; channels whose
% buffer power exceeds a threshold for long enough, and which form a bipolar
% (adjacent-contact) pair, are flagged as the stim pair.

%% Defaults (ported from find_ad_fcn)
d.chunkDuration            = 0.02;   % internal analysis step (s)
d.decay                    = 0.3;    % buffer decay (emphasizes 50 Hz repetition)
d.stimPowerBoost           = 1e9;    % buffer power threshold to call stim on
d.stimOffPower             = 1e7;    % buffer power below this => stim off
d.secs_thresh_stim         = 0.3;    % window over which to require power
d.perc_above_thresh_stim   = 0.5;    % fraction of that window above threshold
d.min_event_sep            = 0.5;    % (s) ignore re-detections within this gap
if nargin < 2 || isempty(params), params = struct(); end
params = set_defaults(params, d);

fs        = session.fs;
chLabels  = session.chLabels;
values    = session.values;
start_t   = session.start_time;
numCh     = numel(chLabels);

chunkSize       = round(params.chunkDuration*fs);
num_thresh_stim = ceil(params.secs_thresh_stim/params.chunkDuration);

% bipolar structure (static for the session)
[~,~,altBipolarIndices] = find_bipolar_pairs(chLabels,1:numCh);
exclude = find_exclude_chs(chLabels);

%% State
buffer         = zeros(chunkSize,numCh);
last_ones_stim = zeros(num_thresh_stim,numCh);
stim_on        = 0;
keep_pair      = [];
last_on_idx    = -inf;

stimT = table('Size',[0 5], ...
    'VariableTypes',{'double','cell','double','double','double'}, ...
    'VariableNames',{'EventNum','StimChs','StimIdx','OnTime','OffTime'});
stimT.StimIdx = zeros(0,2);   % StimIdx holds a [ch1 ch2] pair per row
ev = 0;

numSamples = size(values,1);
last_on_rownum = nan;

for startIdx = 1:chunkSize:(numSamples - chunkSize)
    endIdx = startIdx + chunkSize - 1;
    chunk = values(startIdx:endIdx,:);

    % nan handling + demean
    if all(isnan(chunk),'all'), chunk = zeros(size(chunk)); end
    for ich = 1:numCh
        chunk(isnan(chunk(:,ich)),ich) = mean(chunk(:,ich),'omitnan');
    end
    chunk = chunk - mean(chunk,1,'omitnan');
    chunk(:,exclude) = nan;

    % decaying buffer => repetitive (50 Hz) signal accumulates
    buffer(isnan(buffer)) = chunk(isnan(buffer));
    buffer = buffer*params.decay + chunk;
    buffer_power = sum(buffer.^2,1);

    abs_time = endIdx/fs + start_t;

    if stim_on == 0
        % looking for stim onset
        chs_above = buffer_power > params.stimPowerBoost;
        last_ones_stim(1:end-1,:) = last_ones_stim(2:end,:);
        last_ones_stim(end,:) = chs_above;
        detected = sum(last_ones_stim==1,1) > size(last_ones_stim,1)*params.perc_above_thresh_stim;

        [~,bip_idx] = find_bipolar_pairs(chLabels(detected),find(detected));
        if ~isempty(bip_idx) && (endIdx - last_on_idx) > params.min_event_sep*fs
            % choose highest-power bipolar pair
            if size(bip_idx,1) > 1
                mp = arrayfun(@(j) mean(buffer_power(bip_idx(j,:))), 1:size(bip_idx,1));
                [~,h] = max(mp);
                keep_pair = bip_idx(h,:);
            else
                keep_pair = bip_idx;
            end
            ev = ev + 1;
            % Build the row as a one-row table and append. The double-wrapped
            % StimChs {{a,b}} stores the 1x2 label cell intact; keep_pair
            % (1x2) makes StimIdx a proper 2-wide column.
            newRow = table(ev, {{chLabels{keep_pair(1)}, chLabels{keep_pair(2)}}}, ...
                keep_pair, abs_time, NaN, ...
                'VariableNames', {'EventNum','StimChs','StimIdx','OnTime','OffTime'});
            stimT = [stimT; newRow]; %#ok<AGROW>
            last_on_rownum = height(stimT);
            last_on_idx = endIdx;
            stim_on = 1;
            last_ones_stim(:) = 0;
        end
    else
        % looking for stim offset on the stim pair
        if mean(buffer_power(keep_pair)) < params.stimOffPower
            stim_on = 0;
            stimT.OffTime(last_on_rownum) = abs_time;
        end
    end
end

% close any dangling event at session end
open = isnan(stimT.OffTime);
stimT.OffTime(open) = session.end_time;

end


function s = set_defaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s,f{i}) || isempty(s.(f{i})), s.(f{i}) = d.(f{i}); end
end
end
