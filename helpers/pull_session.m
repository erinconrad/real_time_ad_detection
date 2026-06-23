function session = pull_session(ieeg_name, start_time, end_time, modifier)
% PULL_SESSION  Download a (possibly long) session window from ieeg.org in
% time chunks and return a single session struct.
%
%   session = pull_session(ieeg_name, start_time, end_time, modifier)
%
% ieeg.org rejects any single request where (hz * channels * seconds) exceeds
% 500*130*2000 = 1.3e8. Whole stim sessions can be thousands of seconds long,
% so we split the time range into chunks small enough to stay under that cap
% (with a safety margin) and vertically concatenate the samples.
%
% Returns a struct with fields matching export_sessions:
%   .values .fs .chLabels .start_time .end_time .ieeg_name .modifier .aT

if nargin < 4, modifier = 1; end

L = rt_paths;   % sets up paths (project + IEEG toolbox)
pwfile = L.ieeg_pw_file; login_name = L.ieeg_login;

% --- one tiny pull to get fs, channel labels, and annotations ---
meta = download_ieeg_data(ieeg_name, login_name, pwfile, [start_time start_time+1], 1);
fs   = meta.fs;
nchs = size(meta.chLabels,1);

% --- size the time chunks ---
budget = 500*130*2000;     % ieeg.org cap on hz*channels*seconds
safety = 0.8;              % stay comfortably under the cap
nsamp_chunk = max(round(fs), floor(safety*budget/nchs));   % samples per chunk

samp0 = round(start_time*fs);
samp1 = round(end_time*fs);
nTotal = samp1 - samp0 + 1;
nChunks = ceil(nTotal/nsamp_chunk);

fprintf('Pulling %s_%d in %d chunk(s) (%.0f s, %d ch, %g Hz)...\n', ...
    ieeg_name, modifier, nChunks, end_time-start_time, nchs, fs);

vals = [];
s = samp0; k = 0;
while s <= samp1
    k = k + 1;
    e = min(samp1, s + nsamp_chunk - 1);
    fprintf('  chunk %d/%d: %.1f-%.1f s\n', k, nChunks, s/fs, e/fs);
    d = download_ieeg_data(ieeg_name, login_name, pwfile, [s/fs, e/fs], 0);
    vals = [vals; d.values]; %#ok<AGROW>
    s = e + 1;
end

session = struct();
session.values     = vals;
session.fs         = fs;
session.chLabels   = decompose_labels(meta.chLabels(:,1));
session.start_time = start_time;
session.end_time   = end_time;
session.ieeg_name  = ieeg_name;
session.modifier   = modifier;
if isfield(meta,'aT'), session.aT = meta.aT; else, session.aT = table(); end

end
