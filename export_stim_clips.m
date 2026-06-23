function clipTable = export_stim_clips(session_rows, pre_s, post_s, overwrite, chmode, params, max_clips)
% EXPORT_STIM_CLIPS  Save one clip per individual stimulation. This is the
% main data step (replaces full-session export): clips ARE the working data.
%
%   export_stim_clips                       % all sessions, 10 s pre / 10 s post
%   export_stim_clips(1)                    % only session row 1
%   export_stim_clips([], pre_s, post_s)    % custom window
%   export_stim_clips(rows, pre, post, overwrite, chmode, params, max_clips)
%
%   chmode    : 'electrode' (default) saves every contact on the stim electrode;
%               'neighbors' saves only the stim pair +/- 1 contact.
%   max_clips : stop after saving this many NEW clips (default Inf). Handy for a
%               quick test, e.g. export_stim_clips(1,[],[],[],[],[],2). See also
%               test_export_clips for a one-call test that also plots.
%
% Memory-safe by design: it never loads a whole session. For each session it
%   (1) streams the data in chunks to find the stim events (find_session_stims),
%   (2) for each stim, pulls ONLY the short clip window [stimOn-pre, stimOff+post]
%       and keeps ONLY the stim-electrode channels.
%
% Outputs (clips/):
%   <ieeg>_<mod>_ev<NN>_<LBa-LBb>.mat   struct `clip` (lossless, true times)
%   <ieeg>_<mod>_ev<NN>_<LBa-LBb>.edf   EDF (NaN->0, zero-padded to whole s)
%   clip_index.csv                      master index of all clips
%
% RUN test_stim_detection FIRST to confirm stim detection looks right.

if nargin < 1, session_rows = []; end
if nargin < 2 || isempty(pre_s),  pre_s  = 10; end
if nargin < 3 || isempty(post_s), post_s = 10; end
if nargin < 4 || isempty(overwrite), overwrite = 0; end
if nargin < 5 || isempty(chmode), chmode = 'electrode'; end
if nargin < 6, params = struct(); end
if nargin < 7 || isempty(max_clips), max_clips = Inf; end

L = rt_paths;   % sets up paths (project + IEEG toolbox)
pwfile = L.ieeg_pw_file; login_name = L.ieeg_login;
clip_dir = L.clip_dir;

fT = readtable(L.session_csv);
if isempty(session_rows), session_rows = 1:height(fT); end
vn = fT.Properties.VariableNames;
endcol = vn{find(ismember(lower(vn), {'end','xend'}), 1)};

clipTable = table();
n_saved = 0;

for si = session_rows(:)'
    if n_saved >= max_clips, break; end
    ieeg_name  = fT.ieeg_name{si};
    modifier   = fT.Modifier(si);
    start_time = fT.start(si);
    end_time   = fT.(endcol)(si);
    tag = sprintf('%s_%d', ieeg_name, modifier);

    % --- pass 1: find stims (streamed, nothing kept in memory) ---
    [stimT, fs, chLabels] = find_session_stims(ieeg_name, start_time, end_time, ...
        modifier, [], params);
    if isempty(stimT)
        fprintf('%s: no stim events detected.\n', tag);
        continue;
    end

    % --- pass 2: pull + save a short clip per stim ---
    for e = 1:height(stimT)
        if n_saved >= max_clips, break; end
        stimOn  = stimT.OnTime(e);
        stimOff = stimT.OffTime(e);
        pairStr = strjoin(sort(stimT.StimChs{e}),'-');

        [idx, labs] = stim_clip_channels(chLabels, stimT.StimChs{e}, chmode);
        if isempty(idx)
            warning('%s event %d: no clip channels for %s, skipping.', tag, e, pairStr);
            continue;
        end

        base = regexprep(sprintf('%s_ev%02d_%s', tag, e, pairStr), '[^\w\-]', '');
        mat_path = fullfile(clip_dir, [base '.mat']);
        edf_path = fullfile(clip_dir, [base '.edf']);
        if ~overwrite && exist(mat_path,'file')
            continue;
        end

        % short targeted pull of the clip window (all channels, then subset)
        w0 = stimOn - pre_s; w1 = stimOff + post_s;
        d = download_ieeg_data(ieeg_name, login_name, pwfile, [w0 w1], 0);
        clipVals = d.values(:, idx);

        clip = struct();
        clip.values     = clipVals;
        clip.fs         = fs;
        clip.labels     = labs;
        clip.chmode     = chmode;
        clip.ieeg_name  = ieeg_name;
        clip.modifier   = modifier;
        clip.event      = e;
        clip.stim_pair  = stimT.StimChs{e};
        clip.stimOn     = stimOn;
        clip.stimOff    = stimOff;
        clip.clip_start = w0;
        clip.clip_end   = w1;
        clip.pre_s      = pre_s;
        clip.post_s     = post_s;
        save(mat_path, 'clip');

        try
            write_edf_clip(edf_path, clipVals, fs, labs);
        catch ME
            fprintf('  [EDF skipped for %s: %s]\n', base, ME.message);
            edf_path = '';
        end

        clipTable = [clipTable; table({base}, {tag}, e, {pairStr}, {strjoin(labs,',')}, ...
            stimOn, stimOff, w0, w1, {mat_path}, {edf_path}, ...
            'VariableNames', {'clip_name','session','event','stim_pair','channels', ...
            'stimOn','stimOff','clip_start','clip_end','mat','edf'})]; %#ok<AGROW>

        n_saved = n_saved + 1;
        fprintf('  saved %s  (%d ch, %.1f s)\n', base, numel(idx), w1-w0);
    end
end

if ~isempty(clipTable)
    idx_csv = fullfile(clip_dir, 'clip_index.csv');
    if exist(idx_csv,'file') && ~overwrite
        old = readtable(idx_csv, 'TextType','char');
        clipTable = [old; clipTable];
        clipTable = unique(clipTable, 'rows', 'stable');
    end
    writetable(clipTable, idx_csv);
    fprintf('\nSaved/updated index -> %s (%d clips total)\n', idx_csv, height(clipTable));
end

end
