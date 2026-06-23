function [stimT, cmp] = test_stim_detection(session_row, max_stims, chunk_minutes, params)
% TEST_STIM_DETECTION  Sanity-check automatic stim channel-pair detection on
% ONE session before committing to saving clips for everything.
%
%   [stimT, cmp] = test_stim_detection                 % first session, first 3 stims
%   [stimT, cmp] = test_stim_detection(1, 3)           % session row 1, first 3 stims
%   [stimT, cmp] = test_stim_detection(row, max_stims, chunk_minutes, params)
%
% To stay fast, it does NOT download the whole session. It streams the data
% in `chunk_minutes`-long chunks (default 5 min), runs find_stim_events on
% each chunk, and STOPS as soon as it has found `max_stims` stimulations
% (default 3). Each detected bipolar stim pair is compared against the
% "Closed relay to X and Y" ieeg annotations (the true stim contacts), and
% each found stim is plotted (pair + neighbors, with stim on/off + relay).
%
% Use this to confirm detection quality, then run export_stim_clips.
%
% Note: a stim that straddles a chunk boundary could be split/missed; for a
% quick first-few-stims check this is rarely an issue. Increase chunk_minutes
% if a session starts with a long quiet period before the first stim.

if nargin < 1 || isempty(session_row),   session_row = 1;   end
if nargin < 2 || isempty(max_stims),     max_stims = 3;     end
if nargin < 3 || isempty(chunk_minutes), chunk_minutes = 5; end
if nargin < 4,                           params = struct(); end

L = rt_paths;   % sets up paths (project + IEEG toolbox)
pwfile = L.ieeg_pw_file; login_name = L.ieeg_login;

%% Which session
fT = readtable(L.session_csv);
ieeg_name  = fT.ieeg_name{session_row};
modifier   = fT.Modifier(session_row);
start_time = fT.start(session_row);
vn = fT.Properties.VariableNames;
endcol = vn{find(ismember(lower(vn), {'end','xend'}), 1)};
end_time = fT.(endcol)(session_row);

fprintf('Testing stim detection on %s_%d  (%.1f - %.1f s)\n', ...
    ieeg_name, modifier, start_time, end_time);

%% One tiny pull to get fs, channel labels, and relay annotations
meta = download_ieeg_data(ieeg_name, login_name, pwfile, [start_time start_time+1], 1);
fs       = meta.fs;
chLabels = decompose_labels(meta.chLabels(:,1));
relay    = parse_relay_annotations(meta.aT);

%% Stream chunks until we have the first few stims
chunk_s = chunk_minutes*60;
stimT = []; cmp = table();
t = start_time;
fprintf('Scanning in %g-min chunks for the first %d stim(s)...\n', chunk_minutes, max_stims);
while t < end_time && height(cmp) < max_stims
    t2 = min(end_time, t + chunk_s);
    fprintf('  chunk %.1f-%.1f s ...\n', t, t2);
    d = download_ieeg_data(ieeg_name, login_name, pwfile, [t t2], 0);

    chunkSession = struct('values',d.values,'fs',fs,'chLabels',{chLabels}, ...
        'start_time',t,'end_time',t2,'ieeg_name',ieeg_name,'modifier',modifier,'aT',table());

    sT = find_stim_events(chunkSession, params);
    for e = 1:height(sT)
        dp = sort(sT.StimChs{e});
        [relayStr, relayTime, matchStr] = match_relay(relay, sT.OnTime(e), dp);

        plot_one_stim(chunkSession, sT, e, relayStr);

        if isempty(stimT), stimT = sT(e,:); else, stimT = [stimT; sT(e,:)]; end %#ok<AGROW>
        cmp = [cmp; table(height(cmp)+1, {strjoin(dp,'-')}, sT.OnTime(e), sT.OffTime(e), ...
            sT.OffTime(e)-sT.OnTime(e), {relayStr}, relayTime, string(matchStr), ...
            'VariableNames',{'Event','DetectedPair','OnTime','OffTime','Dur_s', ...
                             'RelayPair','RelayTime','Match'})]; %#ok<AGROW>

        if height(cmp) >= max_stims, break; end
    end
    t = t2;
end

if isempty(stimT)
    warning(['No stim detected in scanned range -- consider lowering ' ...
        'stimPowerBoost / perc_above_thresh_stim, or raising chunk_minutes.']);
    return;
end
stimT.EventNum = (1:height(stimT))';

%% Report
fprintf('\n');
disp(cmp);
scored = ismember(cmp.Match, ["yes","NO"]);
if any(scored)
    fprintf('Stim-pair match vs relay annotations: %d / %d (%.0f%%)\n', ...
        sum(cmp.Match=="yes"), sum(scored), 100*mean(cmp.Match(scored)=="yes"));
else
    fprintf('No relay annotations available to validate against.\n');
end
fprintf('(Scanned %.1f min of data to find %d stim(s).)\n', (t-start_time)/60, height(cmp));

end


% ------------------------------------------------------------------------
function [relayStr, relayTime, matchStr] = match_relay(relay, onTime, dp)
relayStr = ''; relayTime = NaN; matchStr = "no-ann";
if isempty(relay), return; end
dt = relay.time - onTime;             % relay typically closes just before stim
cand = find(dt > -10 & dt < 3);
if isempty(cand), [~,cand] = min(abs(dt)); end
[~,best] = min(abs(dt(cand)));
ri = cand(best);
rp = sort({relay.chA{ri}, relay.chB{ri}});
relayStr = strjoin(rp,'-');
relayTime = relay.time(ri);
if isequal(dp(:), rp(:)), matchStr = "yes"; else, matchStr = "NO"; end
end


% ------------------------------------------------------------------------
function relay = parse_relay_annotations(aT)
relay = struct('time',[],'chA',{{}},'chB',{{}});
if isempty(aT) || ~ismember('Type', aT.Properties.VariableNames)
    relay = []; return;
end
times = []; A = {}; B = {};
for i = 1:height(aT)
    typ = aT.Type{i};
    if ischar(typ) || isstring(typ)
        tok = regexp(char(typ), 'Closed relay to ([A-Za-z]+\d+)\s+and\s+([A-Za-z]+\d+)', 'tokens', 'once');
        if ~isempty(tok)
            times(end+1,1) = aT.Start(i);                 %#ok<AGROW>
            a = decompose_labels({tok{1}}); b = decompose_labels({tok{2}});
            A{end+1,1} = a{1};                            %#ok<AGROW>
            B{end+1,1} = b{1};                            %#ok<AGROW>
        end
    end
end
if isempty(times), relay = []; return; end
relay.time = times; relay.chA = A; relay.chB = B;
end


% ------------------------------------------------------------------------
function plot_one_stim(session, stimT, e, relayStr)
fs = session.fs;
chLabels = session.chLabels;
stimOn  = stimT.OnTime(e);
stimOff = stimT.OffTime(e);

[idx, labs] = stim_clip_channels(chLabels, stimT.StimChs{e});
if isempty(idx), idx = stimT.StimIdx(e,:)'; labs = chLabels(idx); end

pre = 5; post = 8;
i0 = max(1, round((stimOn - pre  - session.start_time)*fs));
i1 = min(size(session.values,1), round((stimOff + post - session.start_time)*fs));
tt = session.start_time + (i0:i1)/fs;
V  = session.values(i0:i1, idx);
V  = V - mean(V,1,'omitnan');

figure('Name',sprintf('Stim event %d  (%s)', e, strjoin(sort(stimT.StimChs{e}),'-')), ...
    'Position',[60 60 1300 700]);
hold on;
sd = nanstd(V(:)); if sd==0 || isnan(sd), sd = 1; end
step = 8*sd;
for c = 1:numel(idx)
    plot(tt, V(:,c) - (c-1)*step, 'k');
    text(tt(1), -(c-1)*step, [labs{c} ' '], 'HorizontalAlignment','right','FontSize',11);
end
yl = [-numel(idx)*step, step]; ylim(yl);
plot([stimOn stimOn], yl, 'b', 'LineWidth',1.5);
plot([stimOff stimOff], yl, 'r', 'LineWidth',1.5);
text(stimOn, yl(2), ' stim on (det)', 'Color','b');
text(stimOff, yl(2), ' stim off (det)', 'Color','r');
ax = gca; ax.XAxis.Exponent = 0; set(ax,'YTick',[]);
xlabel('Time (s)');
title(sprintf('Event %d: detected %s | relay annotation: %s', ...
    e, strjoin(sort(stimT.StimChs{e}),'-'), relayStr), 'Interpreter','none');
end
