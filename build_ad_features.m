function build_ad_features(params, remote_clips)
% BUILD_AD_FEATURES  Pre-compute and cache detector features for every
% annotated clip, so threshold/CV tuning in validate_ad_detector is fast and
% does not need repeated clip access.
%
%   build_ad_features                       % local clips/ dir
%   build_ad_features(params)               % override detector params
%   build_ad_features(params, 'user@host:/path/to/clips')   % stream clips
%
% Reads ground_truth/clip_gt.csv, computes ad_clip_features for each clip with
% AD in {y,n}, and saves results/ad_features.mat (features + labels + grouping
% + the annotation fields needed for the exclusion rules).

if nargin < 1, params = struct(); end
if nargin < 2, remote_clips = ''; end
L = rt_paths;            % sets up the path (incl. helpers/) first
p = ad_params(params);

gt_csv = fullfile(L.gt_dir, 'clip_gt.csv');
assert(exist(gt_csv,'file')==2, 'No ground truth at %s', gt_csv);
gt = readtable(gt_csv, 'TextType','char');
AD = string(gt.AD);
sel = find(ismember(AD, ["y","n"]));
assert(~isempty(sel), 'No annotated (y/n) clips in %s', gt_csv);

remote = ~isempty(remote_clips);
if remote
    cache = fullfile(tempdir,'ad_clip_cache'); if ~exist(cache,'dir'), mkdir(cache); end
end

feats = struct('clip_name',{},'patient',{},'label',{},'AD',{}, ...
    'AD_onset',{},'AD_offset',{},'ad_ongoing',{},'stimOn',{},'stimOff',{}, ...
    'clip_start',{},'clip_end',{},'F',{});

% session -> sorted stim onset times (for the next-stim cap)
stimMap = load_stim_times(L, remote_clips);

fprintf('Building features for %d annotated clips...\n', numel(sel));
for k = 1:numel(sel)
    i = sel(k);
    cname = char(gt.clip_name(i));

    if remote
        local_mat = fullfile(cache,[cname '.mat']);
        if ~scp_fetch(remote_clips,[cname '.mat'],local_mat)
            warning('fetch failed %s', cname); continue;
        end
    else
        local_mat = fullfile(L.clip_dir,[cname '.mat']);
        if ~exist(local_mat,'file')
            warning('missing clip %s', cname); continue;
        end
    end
    C = load(local_mat); clip = C.clip;

    % next stim onset in this session (from the stim event list)
    tag = sprintf('%s_%d', clip.ieeg_name, clip.modifier);
    if isKey(stimMap, tag)
        ons = stimMap(tag);
        nxt = ons(ons > clip.stimOff + 1e-6);
        if ~isempty(nxt), clip.next_stim_on = nxt(1); else, clip.next_stim_on = inf; end
    elseif ~isfield(clip,'next_stim_on')
        clip.next_stim_on = inf;
    end

    F = ad_clip_features(clip, p);

    r.clip_name  = cname;
    r.patient    = patient_of(clip.ieeg_name);
    r.label      = strcmp(char(gt.AD(i)),'y');
    r.AD         = char(gt.AD(i));
    r.AD_onset   = todbl(gt.AD_onset(i));
    r.AD_offset  = todbl(gt.AD_offset(i));
    r.ad_ongoing = todbl(getcol(gt,'ad_ongoing',i));
    r.stimOn     = clip.stimOn;
    r.stimOff    = clip.stimOff;
    r.clip_start = clip.clip_start;
    r.clip_end   = clip.clip_end;
    r.F          = F;
    feats(end+1) = r; %#ok<AGROW>

    if remote && exist(local_mat,'file'), delete(local_mat); end
    if mod(k,25)==0, fprintf('  %d/%d\n', k, numel(sel)); end
end

out = fullfile(L.results_dir,'ad_features.mat');
save(out, 'feats', 'p', '-v7.3');
fprintf('Saved %d feature sets -> %s\n', numel(feats), out);
end


function pt = patient_of(ieeg_name)
m = regexp(char(ieeg_name), 'HUP\d+', 'match', 'once');
if isempty(m), pt = char(ieeg_name); else, pt = m; end
end

function d = todbl(x)
if isnumeric(x), d = double(x); else, d = str2double(string(x)); end
end

function v = getcol(t, name, i)
if ismember(name, t.Properties.VariableNames), v = t.(name)(i); else, v = NaN; end
end


function M = load_stim_times(L, remote_clips)
% Map: session tag -> sorted stim onset times for ALL stims in the session.
% Built by scanning the clip .mat files (robust; the clip_index.csv can be
% malformed). Cached to results/stim_times.mat -- delete that file to refresh
% after exporting new clips.
M = containers.Map('KeyType','char','ValueType','any');

cache = fullfile(L.results_dir,'stim_times.mat');
if exist(cache,'file')==2
    S = load(cache); M = S.M; return;
end

if ~isempty(remote_clips)
    % don't scp the whole clip set just for times; skip (cap disabled)
    warning('Remote build: next-stim cap disabled (run build on the server to enable).');
    return;
end

d = dir(fullfile(L.clip_dir,'*.mat'));
if isempty(d)
    warning('No clips found to build stim-time index; next-stim cap disabled.');
    return;
end

fprintf('Building stim-time index by scanning %d clips (one-time)...\n', numel(d));
acc = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(d)
    try
        C = load(fullfile(L.clip_dir, d(i).name), 'clip'); clip = C.clip;
    catch
        continue;
    end
    tag = sprintf('%s_%d', clip.ieeg_name, clip.modifier);
    if isKey(acc, tag), acc(tag) = [acc(tag); clip.stimOn]; else, acc(tag) = clip.stimOn; end
end
k = keys(acc);
for i = 1:numel(k), M(k{i}) = sort(acc(k{i})); end

save(cache, 'M');
fprintf('Saved stim-time index -> %s\n', cache);
end


% ---- minimal scp helpers (mirror annotate_ground_truth) ----
function [host, rpath] = split_remote(remote)
kk = strfind(remote, ':'); kk = kk(1);
host = remote(1:kk-1); rpath = regexprep(remote(kk+1:end), '/+$', '');
end

function ok = scp_fetch(remote, fname, local_path)
[host, rpath] = split_remote(remote);
cmd = sprintf('scp -q "%s:%s/%s" "%s"', host, rpath, fname, local_path);
ok = (system(cmd)==0) && exist(local_path,'file')==2;
end
