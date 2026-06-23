function cfg = rt_local_config
% RT_LOCAL_CONFIG  Per-machine settings for the ad_realtime_prelim project.
%
% SETUP ON A NEW MACHINE / SERVER:
%   1. Copy this file to  rt_local_config.m  (same folder, drop "_template").
%   2. Edit the three paths below for this machine.
%   3. rt_local_config.m is gitignored, so your credentials are never committed.
%
% rt_paths.m uses this file if it exists. If it does not exist, rt_paths falls
% back to seizure_termination_paths (the original lab paths) for backwards
% compatibility on the original machine.
%
% Only these three fields are required by the pipeline:

% Folder containing the IEEG MATLAB toolbox (the one with IEEGSession.m)
cfg.ieeg_folder  = '/path/to/ieeg-matlab-1.13.2';

% Your ieeg.org password .bin file (created with IEEGSession's pwgen tool)
cfg.ieeg_pw_file = '/path/to/your_ieeglogin.bin';

% Your ieeg.org username
cfg.ieeg_login   = 'your_ieeg_username';

end
