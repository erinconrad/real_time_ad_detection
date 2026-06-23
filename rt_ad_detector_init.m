function state = rt_ad_detector_init(chLabels, fs, t0, file_name, params)
% RT_AD_DETECTOR_INIT  Initialize the streaming AD detector state.
%
%   state = rt_ad_detector_init(chLabels, fs, t0, file_name)
%   state = rt_ad_detector_init(chLabels, fs, t0, file_name, params)
%
% This detector is fully causal: feed it data one block at a time with
% rt_ad_detector_step (blocks may be any length, e.g. 100-1000 ms). It
% detects 50 Hz stim (channel pair + on/off) and afterdischarges online,
% maintaining all internal state (decaying buffer, rolling baselines,
% high-pass filter state, threshold counters) across blocks.
%
% Parameters default to the values tuned in the original find_ad_fcn.

if nargin < 4, file_name = ''; end
if nargin < 5, params = struct(); end

% ---- default parameters (from find_ad_fcn) ----
d.chunkDuration          = 0.02;   % internal analysis window (s)
d.updateInterval         = 0.02;   % advance per internal step (s)
d.hpf_alpha              = 0.99;   % stevefilter leak coefficient

% stim detection
d.decay                  = 0.3;
d.stimPowerBoost         = 1e9;
d.stimOffPower           = 1e7;
d.secs_thresh_stim       = 0.3;
d.perc_above_thresh_stim = 0.5;

% AD detection
d.n_baseline_all         = 100;
d.n_baseline_keep        = 50;
d.stopLookingADSecs      = 5;
d.ad_thresh              = 30;
d.ad_too_high_thresh     = 1e4;
d.coolDownLastSat        = 2;
d.secs_thresh            = 2;
d.num_above_thresh       = 10;
d.hfband                 = [400 500];
d.hfthresh               = 1e4;

% bad channels
d.bad_ch_amp             = 1e4;
d.n_bad                  = 10;
d.n_bad_reset            = 100;
d.n_reduce_bad           = 5;

p = set_defaults(params, d);

numCh = numel(chLabels);

state.p          = p;
state.fs         = fs;
state.t0         = t0;          % absolute time (s) of sample index 0
state.chLabels   = chLabels(:);
state.numChannels= numCh;
state.file_name  = file_name;

% derived sizes
state.chunkSize       = round(p.chunkDuration*fs);
state.updateSize      = round(p.updateInterval*fs);
state.num_thresh      = round(p.secs_thresh/p.chunkDuration);
state.num_thresh_stim = ceil(p.secs_thresh_stim/p.chunkDuration);

% static channel structure
[~,~,state.altBipolarIndices] = find_bipolar_pairs(state.chLabels,1:numCh);
state.exclude = find_exclude_chs(state.chLabels);

% streaming buffers
state.residual    = zeros(0,numCh); % unprocessed tail samples
state.samp_offset = 0;              % global index (0-based) of residual(1)
state.zi          = zeros(1,numCh); % HPF state

% detector state
state.buffer         = zeros(state.chunkSize,numCh);
state.last_ones      = zeros(state.num_thresh,numCh);
state.last_ones_stim = zeros(state.num_thresh_stim,numCh);
state.baselines_all  = repmat({zeros(state.chunkSize,numCh)}, p.n_baseline_all, 1);
state.baselines      = state.baselines_all(1:p.n_baseline_keep);

state.bad_ch_counter      = zeros(1,numCh);
state.n_bad_loop_counter  = 0;
state.time_since_last_sat = nan(1,numCh);

state.stim_on        = 0;
state.look_for_sat   = 0;
state.last_stim_chs  = [];
state.keep_pair      = [];
state.last_stim_on   = NaN;
state.last_stim_off  = inf;
state.last_stim_rows = [];
state.ad_chs_this_stim = [];

% master event log
state.events = empty_event_table();

end


function T = empty_event_table()
T = table('Size',[0 4], ...
    'VariableTypes',{'cell','cell','double','double'}, ...
    'VariableNames',{'Type','Channels','OnTime','OffTime'});
end

function s = set_defaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s,f{i}) || isempty(s.(f{i})), s.(f{i}) = d.(f{i}); end
end
end
