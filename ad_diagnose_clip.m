function ad_diagnose_clip(clip_name, params, remote_clips)
% AD_DIAGNOSE_CLIP  Explain why a clip was / wasn't detected. Plots the
% candidate-channel line-length z-scores per window with the threshold, and
% marks which windows were vetoed and why (HF guard, 60 Hz line exclusion, or
% next-stim cap). Useful for chasing false negatives / false positives.
%
%   ad_diagnose_clip('HUP218_phaseII_D02_1_ev02_RI4-RI5')
%   ad_diagnose_clip(name, struct('T',1.62))
%   ad_diagnose_clip(name, p, 'user@host:/path/to/clips')

if nargin < 2, params = struct(); end
if nargin < 3, remote_clips = ''; end
L = rt_paths;
p = ad_params(params);

% load clip
if ~isempty(remote_clips)
    lp = fullfile(tempdir,[clip_name '.mat']);
    [host,rpath] = split_remote(remote_clips);
    system(sprintf('scp -q "%s:%s/%s.mat" "%s"', host, rpath, clip_name, lp));
else
    lp = fullfile(L.clip_dir,[clip_name '.mat']);
end
assert(exist(lp,'file')==2, 'Clip not found: %s', lp);
C = load(lp); clip = C.clip;

% next-stim time from clip_index (so the cap matches validation)
if ~isfield(clip,'next_stim_on') || isempty(clip.next_stim_on)
    clip.next_stim_on = next_stim_from_index(L, clip);
end

[F, dbg] = ad_clip_features(clip, p);
[score, onset, chan] = ad_apply_rule(F, p.N);
detected = score > p.T;

% --- report ---
fprintf('\n%s\n', clip_name);
fprintf('candidate channels: %s\n', strjoin(F.chan_labels, ', '));
fprintf('threshold T=%.2f, need N=%d consecutive valid windows\n', p.T, p.N);
fprintf('best run min-z = %.2f on %s  => %s\n', score, chan, ...
    ternary(detected,'DETECTED','not detected'));
for j = 1:numel(F.chan_labels)
    zc = F.z(:,j);
    fprintf('  %-10s maxz=%.1f  HFveto=%d  lineExcl=%d\n', F.chan_labels{j}, ...
        max(zc), sum(dbg.hf_bad(:,j)), dbg.line_bad(j));
end
if any(dbg.cap_bad)
    fprintf('  next-stim cap removed %d of %d windows (from %.1fs)\n', ...
        sum(dbg.cap_bad), numel(dbg.cap_bad), F.win_t(find(dbg.cap_bad,1)));
end

% --- plot z per channel ---
figure('Name',['diagnose ' clip_name],'Position',[60 60 1200 700],'Color','w');
hold on;
cols = lines(numel(F.chan_labels));
for j = 1:numel(F.chan_labels)
    valid = ~dbg.hf_bad(:,j) & ~dbg.line_bad(j) & ~dbg.cap_bad;
    plot(F.win_t, F.z(:,j), '-o', 'Color', cols(j,:), 'MarkerFaceColor', cols(j,:), ...
        'DisplayName', F.chan_labels{j});
    % mark vetoed windows with an x
    bad = ~valid;
    if any(bad)
        plot(F.win_t(bad), F.z(bad,j), 'x', 'Color',[0 0 0], 'MarkerSize',10, ...
            'LineWidth',1.5, 'HandleVisibility','off');
    end
end
yl = ylim;
plot(xlim, [p.T p.T], 'r--', 'DisplayName', sprintf('T=%.2f',p.T));
% shade next-stim-capped region
if any(dbg.cap_bad)
    ct = F.win_t(find(dbg.cap_bad,1));
    patch([ct max(F.win_t)+p.win_s max(F.win_t)+p.win_s ct], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.6 0.6 0.6], 'FaceAlpha',0.12,'EdgeColor','none','DisplayName','next-stim cap');
end
xlabel('Time (s)'); ylabel('line-length z vs baseline');
title(sprintf('%s   (black x = vetoed window)', clip_name),'Interpreter','none');
legend('Location','northeast'); grid on;
end


function nso = next_stim_from_index(L, clip)
nso = inf;
idx = fullfile(L.clip_dir,'clip_index.csv');
if exist(idx,'file')~=2, return; end
T = readtable(idx,'TextType','char'); vn = T.Properties.VariableNames;
if ismember('session',vn), sess = string(T.session);
elseif ismember('clip_name',vn), sess = regexprep(string(T.clip_name),'_ev\d+.*$','');
else, return; end
if ~ismember('stimOn',vn), return; end
on = T.stimOn; if iscell(on), on = str2double(string(on)); end
tag = sprintf('%s_%d', clip.ieeg_name, clip.modifier);
ons = on(sess==tag);
nxt = ons(ons > clip.stimOff + 1e-6);
if ~isempty(nxt), nso = min(nxt); end
end

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end

function [host,rpath]=split_remote(remote)
k=strfind(remote,':'); k=k(1); host=remote(1:k-1); rpath=regexprep(remote(k+1:end),'/+$','');
end
