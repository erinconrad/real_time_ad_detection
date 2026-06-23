function cfg = rt_local_config
% RT_LOCAL_CONFIG  Per-machine settings for the ad_realtime_prelim project.
%
% SETUP ON A NEW MACHINE / SERVER:
%   1. Copy this file to  rt_local_config.m  and edit the three paths below.
%   2. Put rt_local_config.m wherever you keep machine-local config -- it does
%      NOT have to live in the repo. The recommended pattern is to keep it in
%      your own tools/ directory (outside any git repo) and addpath that dir at
%      MATLAB startup. rt_paths finds rt_local_config anywhere on the path.
%      (If you do keep it inside scripts/, it is gitignored.)
%
% rt_paths.m uses rt_local_config if it is on the path. If it is not found,
% rt_paths falls back to seizure_termination_paths (the original lab config),
% so you can also just keep using that on each server with these fields set.
%
% Only these three fields are required by the pipeline:

% Folder containing the IEEG MATLAB toolbox (the one with IEEGSession.m)
cfg.ieeg_folder  = '/path/to/ieeg-matlab-1.13.2';

% Your ieeg.org password .bin file (created with IEEGSession's pwgen tool)
cfg.ieeg_pw_file = '/path/to/your_ieeglogin.bin';

% Your ieeg.org username
cfg.ieeg_login   = 'your_ieeg_username';

end
