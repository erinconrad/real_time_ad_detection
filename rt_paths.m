function L = rt_paths
% RT_PATHS  Paths/config for the real-time AD detection prelim project.
%
% Self-contained and portable: it locates itself relative to this file, so the
% scripts/ folder can be cloned to another machine and run as-is. The ONLY
% per-machine settings (IEEG toolbox folder, ieeg.org password file, login)
% come from rt_local_config.m, which you create from rt_local_config_template.m
% on each machine (and which is gitignored).
%
% The git repo is scripts/ (code + hfs_sessions.csv + README). Generated data
% lives in the PARENT of scripts/, i.e. outside the repo:
%   <workspace>/
%     scripts/        <- the git repo (this file lives here, + hfs_sessions.csv)
%     clips/          <- per-stim clips (the working data)   [outside repo]
%     results/        <- detections, evaluation output       [outside repo]
%     ground_truth/   <- human AD annotations                [outside repo]
%     data/           <- optional full-session .mat          [outside repo]

% ---- per-machine config (toolbox + ieeg credentials) ----
if exist('rt_local_config','file') == 2
    cfg = rt_local_config;                 % preferred: project-local config
elseif exist('seizure_termination_paths','file') == 2
    cfg = seizure_termination_paths;       % fallback: original lab paths
else
    error(['No machine config found. Copy rt_local_config_template.m to ' ...
           'rt_local_config.m (in scripts/) and set ieeg_folder, ' ...
           'ieeg_pw_file, and ieeg_login for this machine.']);
end

L = struct();
L.ieeg_folder  = cfg.ieeg_folder;
L.ieeg_pw_file = cfg.ieeg_pw_file;
L.ieeg_login   = cfg.ieeg_login;

% Code root = folder containing this file (.../ad_realtime_prelim/scripts)
L.rt_root = fileparts(mfilename('fullpath'));

% Project root = parent of scripts/  (.../ad_realtime_prelim)
L.proj_root = fileparts(L.rt_root);

% Generated data/results/clips live under the project root (not in scripts/)
L.out_root    = L.proj_root;
L.data_dir    = fullfile(L.proj_root, 'data');         % optional full sessions
L.results_dir = fullfile(L.proj_root, 'results');      % detections, eval output
L.gt_dir      = fullfile(L.proj_root, 'ground_truth'); % human annotations
L.clip_dir    = fullfile(L.proj_root, 'clips');        % per-stim clips (working data)

% Input session list (lives inside the scripts repo so it is versioned and
% travels with a clone)
L.session_csv = fullfile(L.rt_root, 'hfs_sessions.csv');

% Make sure folders exist
cellfun(@(d) ~exist(d,'dir') && mkdir(d), ...
    {L.data_dir, L.results_dir, L.gt_dir, L.clip_dir});

% Put project scripts + helpers (incl. bundled download_ieeg_data) and the
% IEEG toolbox on the path, so callers don't need their own addpath calls.
addpath(genpath(L.rt_root));
if ~isempty(L.ieeg_folder) && exist(L.ieeg_folder,'dir')
    addpath(genpath(L.ieeg_folder));
end

end
