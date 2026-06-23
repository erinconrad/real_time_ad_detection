function [stimT, fs, chLabels] = find_session_stims(ieeg_name, start_time, end_time, modifier, chunk_minutes, params, merge_gap_s)
% FIND_SESSION_STIMS  Detect all stim events in a session WITHOUT ever
% loading the whole session into memory.
%
%   [stimT, fs, chLabels] = find_session_stims(ieeg_name, start_time, end_time, modifier)
%   [...] = find_session_stims(..., chunk_minutes, params, merge_gap_s)
%
% It streams the session in `chunk_minutes`-long time chunks (default 5),
% runs find_stim_events on each, then merges events that are the same stim
% pair and close together in time (default within merge_gap_s = 3 s). The
% merge cleans up (a) a stim that straddles a chunk boundary and (b) brief
% re-detections of one stimulation.
%
% OUTPUT
%   stimT    : table EventNum, StimChs {1x2}, StimIdx [1x2], OnTime, OffTime (abs s)
%   fs       : sampling rate
%   chLabels : cleaned channel labels (decompose_labels), in dataset order

if nargin < 4 || isempty(modifier),      modifier = 1;      end
if nargin < 5 || isempty(chunk_minutes), chunk_minutes = 5; end
if nargin < 6,                           params = struct(); end
if nargin < 7 || isempty(merge_gap_s),   merge_gap_s = 3;   end

L = rt_paths;   % sets up paths (project + IEEG toolbox)
pwfile = L.ieeg_pw_file; login_name = L.ieeg_login;

% tiny pull for fs + labels
meta = download_ieeg_data(ieeg_name, login_name, pwfile, [start_time start_time+1], 1);
fs       = meta.fs;
chLabels = decompose_labels(meta.chLabels(:,1));

chunk_s = chunk_minutes*60;
raw = empty_stim_table();

t = start_time;
fprintf('  scanning %s_%d for stims in %g-min chunks...\n', ieeg_name, modifier, chunk_minutes);
while t < end_time
    t2 = min(end_time, t + chunk_s);
    d = download_ieeg_data(ieeg_name, login_name, pwfile, [t t2], 0);
    chunkSession = struct('values',d.values,'fs',fs,'chLabels',{chLabels}, ...
        'start_time',t,'end_time',t2,'ieeg_name',ieeg_name,'modifier',modifier,'aT',table());
    sT = find_stim_events(chunkSession, params);
    if ~isempty(sT), raw = [raw; sT]; end %#ok<AGROW>
    t = t2;
end

stimT = merge_stim_events(raw, merge_gap_s);
if ~isempty(stimT)
    stimT.EventNum = (1:height(stimT))';
end
fprintf('  found %d stim event(s).\n', height(stimT));

end


function T = empty_stim_table()
T = table('Size',[0 5], ...
    'VariableTypes',{'double','cell','double','double','double'}, ...
    'VariableNames',{'EventNum','StimChs','StimIdx','OnTime','OffTime'});
T.StimIdx = zeros(0,2);
end


function out = merge_stim_events(T, gap_s)
% Merge rows that are the same stim pair and start within gap_s of the
% previous row's offset (handles boundary splits + re-detections).
out = empty_stim_table();
if isempty(T), return; end

[~,ord] = sort(T.OnTime);
T = T(ord,:);

cur = T(1,:);
for i = 2:height(T)
    same_pair = isequal(sort(T.StimIdx(i,:)), sort(cur.StimIdx(end,:)));
    gap = T.OnTime(i) - cur.OffTime(end);
    if same_pair && gap < gap_s
        % extend current event's offset
        cur.OffTime(end) = max(cur.OffTime(end), T.OffTime(i));
    else
        out = [out; cur]; %#ok<AGROW>
        cur = T(i,:);
    end
end
out = [out; cur];
end
