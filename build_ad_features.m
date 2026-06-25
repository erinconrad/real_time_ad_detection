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
p = ad_params(params);
L = rt_paths;

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
