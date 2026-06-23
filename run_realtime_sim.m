function results = run_realtime_sim(block_ms, params, which_clips)
% RUN_REALTIME_SIM  Replay per-stim CLIPS through the streaming detector one
% block at a time, emulating a real-time data feed.
%
%   results = run_realtime_sim                 % 250 ms blocks, all clips
%   results = run_realtime_sim(block_ms)       % e.g. 100, 250, 500, 1000
%   results = run_realtime_sim(block_ms, params)
%   results = run_realtime_sim(block_ms, params, {'HUP223_phaseII_1_ev01_LB3-LB4'})
%
% For each clip in clips/ it:
%   * loads the clip (.mat, stim-electrode channels only)
%   * feeds the detector blocks of `block_ms` milliseconds
%   * times each block (wall-clock) to verify it keeps up with real time
%   * saves results/<clip>_rt_<block_ms>ms.mat with detected events, the
%     clip's stim times (for evaluation), per-block latency and real_time_factor
%
% The detector auto-detects the stim pair and on/off within the clip (the
% whole electrode is present), then looks for ADs causally afterward.

if nargin < 1 || isempty(block_ms), block_ms = 250;   end
if nargin < 2,                      params = struct(); end
if nargin < 3,                      which_clips = {};  end

L = rt_paths;

listing = dir(fullfile(L.clip_dir, '*.mat'));
names = erase({listing.name}, '.mat');
if ~isempty(which_clips)
    listing = listing(ismember(names, which_clips));
end
assert(~isempty(listing), 'No clips found in %s. Run export_stim_clips first.', L.clip_dir);

results = struct('name',{},'events',{},'stimOn',{},'stimOff',{},'block_ms',{}, ...
    'mean_block_cpu_ms',{},'max_block_cpu_ms',{},'real_time_factor',{});

for i = 1:numel(listing)
    C = load(fullfile(L.clip_dir, listing(i).name));
    clip = C.clip;
    tag = erase(listing(i).name,'.mat');

    fs    = clip.fs;
    vals  = clip.values;
    nSamp = size(vals,1);
    blockSize = max(1, round(block_ms/1000*fs));

    state = rt_ad_detector_init(clip.labels, fs, clip.clip_start, ...
        clip.ieeg_name, params);

    nBlocks = ceil(nSamp/blockSize);
    block_cpu = nan(nBlocks,1);

    b = 0;
    for s0 = 1:blockSize:nSamp
        b = b + 1;
        s1 = min(nSamp, s0 + blockSize - 1);
        tstart = tic;
        [state, ~, ~] = rt_ad_detector_step(state, vals(s0:s1, :)); %#ok<ASGLU>
        block_cpu(b) = toc(tstart)*1000;
    end

    data_seconds = nSamp/fs;
    cpu_seconds  = sum(block_cpu)/1000;

    r.name              = tag;
    r.events            = state.events;
    r.stimOn            = clip.stimOn;
    r.stimOff           = clip.stimOff;
    r.block_ms          = block_ms;
    r.mean_block_cpu_ms = mean(block_cpu,'omitnan');
    r.max_block_cpu_ms  = max(block_cpu);
    r.real_time_factor  = data_seconds / cpu_seconds;
    results(end+1) = r; %#ok<AGROW>

    out = fullfile(L.results_dir, sprintf('%s_rt_%dms.mat', tag, block_ms));
    events = state.events;            %#ok<NASGU>
    meta   = rmfield(r,'events');     %#ok<NASGU>
    save(out, 'events', 'meta');

    n_ad = sum(strcmp(state.events.Type,'AD'));
    fprintf('%s | %d AD detection(s) | block cpu mean %.2f ms, max %.2f ms | %.0fx real time\n', ...
        tag, n_ad, r.mean_block_cpu_ms, r.max_block_cpu_ms, r.real_time_factor);
end

fprintf('\n=== Real-time feasibility (block = %d ms, %d clips) ===\n', block_ms, numel(results));
fprintf('Max single-block CPU: %.2f ms (budget = %d ms)\n', ...
    max([results.max_block_cpu_ms]), block_ms);
fprintf('Median real-time factor: %.0fx faster than real time\n', ...
    median([results.real_time_factor]));

end
