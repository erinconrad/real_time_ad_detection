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
% Map: session tag -> sorted stim onset times, from clip_index.csv (all stims).
M = containers.Map('KeyType','char','ValueType','any');
idx = fullfile(L.clip_dir,'clip_index.csv');
if exist(idx,'file')~=2 && ~isempty(remote_clips)
    tmp = fullfile(tempdir,'clip_index.csv');
    [host,rpath] = split_remote(remote_clips);
    if system(sprintf('scp -q "%s:%s/clip_index.csv" "%s"',host,rpath,tmp))==0
        idx = tmp;
    end
end
if exist(idx,'file')~=2
    warning('No clip_index.csv found; next-stim cap disabled.');
    return;
end
T = readtable(idx,'TextType','char');
vn = T.Properties.VariableNames;

% session key per row (use 'session' column, else derive from clip_name)
if ismember('session', vn)
    sess = string(T.session);
elseif ismember('clip_name', vn)
    sess = regexprep(string(T.clip_name), '_ev\d+.*$', '');
else
    warning('clip_index.csv has no session/clip_name; next-stim cap disabled.');
    M = containers.Map('KeyType','char','ValueType','any'); return;
end

if ~ismember('stimOn', vn)
    warning('clip_index.csv has no stimOn; next-stim cap disabled.');
    M = containers.Map('KeyType','char','ValueType','any'); return;
end
on = T.stimOn;
if iscell(on), on = str2double(string(on)); end

u = unique(sess);
for i = 1:numel(u)
    M(char(u(i))) = sort(on(sess==u(i)));
end
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
