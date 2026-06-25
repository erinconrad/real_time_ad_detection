function review_detections(params, remote_clips, start_filter)
% REVIEW_DETECTIONS  Navigable GUI to review detector output against ground
% truth, one clip at a time, color-coded TP / FP / FN / TN.
%
%   review_detections                         % uses cached features + tuned T
%   review_detections(struct('T',1.62))       % set the decision threshold
%   review_detections(p, 'user@host:/clips')  % stream clips for plotting
%   review_detections(p, '', 'FP')            % start filtered to false positives
%
% Classification uses the cached features (build_ad_features) at threshold T
% (if T not given, the suggested threshold from ad_validation.mat is used).
% The clip signal is loaded only for the clip on screen.
%
% Overlays: stim pair in RED; ground-truth AD channels GREEN (with the AD
% interval as green dashed lines); the detector's onset as a BLUE line on the
% detected channel. The title and its color give the category.
%
% Controls: Prev / Next buttons (or left/right arrows), a category filter
% dropdown, and Quit.

if nargin < 1, params = struct(); end
if nargin < 2, remote_clips = ''; end
if nargin < 3 || isempty(start_filter), start_filter = 'All'; end

L = rt_paths;
p = ad_params(params);

% threshold: explicit > validation-suggested > default
if ~isfield(params,'T')
    vf = fullfile(L.results_dir,'ad_validation.mat');
    if exist(vf,'file'), v = load(vf); if isfield(v,'res'), p.T = v.res.T_suggested; end; end
end

% cached features
ff = fullfile(L.results_dir,'ad_features.mat');
assert(exist(ff,'file')==2, 'Run build_ad_features first (no %s).', ff);
Sl = load(ff); feats = Sl.feats; n = numel(feats);

% gt for AD channels + stim pair string
gt = readtable(fullfile(L.gt_dir,'clip_gt.csv'),'TextType','char');
gt_name = string(gt.clip_name);

% classify every clip
items = struct('clip_name',{},'cat',{},'score',{},'onset',{},'chan',{}, ...
    'gtAD',{},'AD_onset',{},'AD_offset',{},'ad_ongoing',{},'gt_chans',{}, ...
    'stimOn',{},'stimOff',{},'det_t0',{},'det_t1',{},'shorts',{});
for i = 1:n
    [sc, on, ch] = ad_apply_rule(feats(i).F, p.N);
    pred = sc > p.T;
    det_t0 = NaN; det_t1 = NaN;
    if pred
        % contiguous above-threshold span on the detected channel
        F = feats(i).F;
        ci = find(strcmp(F.chan_labels, ch), 1);
        ab = F.z(:,ci) > p.T & F.valid(:,ci);
        k = find(abs(F.win_t - on) < 1e-6, 1);
        if isempty(k), k = find(ab,1); end
        a = k; while a > 1 && ab(a-1), a = a-1; end
        b = k; while b < numel(ab) && ab(b+1), b = b+1; end
        det_t0 = F.win_t(a);
        det_t1 = F.win_t(b) + p.win_s;
    else
        on = NaN; ch = '';   % only report when actually detected
    end

    % supra-threshold runs SHORTER than N windows (above threshold but <3 s):
    % not detections, but worth showing as near-misses
    shorts = struct('t0',{},'t1',{},'chan',{});
    F = feats(i).F;
    if ~isempty(F.win_t)
        for ci = 1:numel(F.chan_labels)
            ab = F.z(:,ci) > p.T & F.valid(:,ci);
            rr = true_runs(ab);
            for q = 1:size(rr,1)
                Lwin = rr(q,2) - rr(q,1) + 1;
                if Lwin >= 1 && Lwin < p.N
                    sh.t0 = F.win_t(rr(q,1));
                    sh.t1 = F.win_t(rr(q,2)) + p.win_s;
                    sh.chan = F.chan_labels{ci};
                    shorts(end+1) = sh; %#ok<AGROW>
                end
            end
        end
    end

    lab  = feats(i).label;

    % exclusions (same rules as validate)
    excl = false;
    if lab
        if feats(i).AD_onset < feats(i).stimOn, excl = true; end
        dur = feats(i).AD_offset - feats(i).AD_onset;
        if ~feats(i).ad_ongoing && dur < p.min_ad_dur, excl = true; end
    end
    if excl
        cat = 'excluded';
    elseif lab && pred,  cat = 'TP';
    elseif ~lab && pred, cat = 'FP';
    elseif lab && ~pred, cat = 'FN';
    else,                cat = 'TN';
    end

    r.clip_name = feats(i).clip_name; r.cat = cat;
    r.score = sc; r.onset = on; r.chan = ch;
    r.gtAD = lab; r.AD_onset = feats(i).AD_onset; r.AD_offset = feats(i).AD_offset;
    r.ad_ongoing = feats(i).ad_ongoing;
    j = find(gt_name==string(feats(i).clip_name),1);
    if ~isempty(j) && ismember('AD_channels',gt.Properties.VariableNames)
        r.gt_chans = char(string(gt.AD_channels(j)));
    else
        r.gt_chans = '';
    end
    r.stimOn = feats(i).stimOn; r.stimOff = feats(i).stimOff;
    r.det_t0 = det_t0; r.det_t1 = det_t1; r.shorts = shorts;
    items(end+1) = r; %#ok<AGROW>
end

% --- figure/state ---
S.items = items; S.p = p; S.clipdir = L.clip_dir;
S.remote = ~isempty(remote_clips); S.remote_clips = remote_clips;
if S.remote
    S.cache = fullfile(tempdir,'ad_review_cache'); if ~exist(S.cache,'dir'), mkdir(S.cache); end
end
S.fig = figure('Name','Review detections','Position',[40 40 1550 860],'Color','w');
S.ax  = axes('Parent',S.fig,'Position',[0.06 0.10 0.78 0.82]);
S.pos = 1; S.gain = 1;

uicontrol(S.fig,'Style','text','Units','normalized','String','Filter:', ...
    'Position',[0.855 0.92 0.14 0.025],'BackgroundColor','w','HorizontalAlignment','left');
S.hFilter = uicontrol(S.fig,'Style','popupmenu','Units','normalized', ...
    'String',{'All','TP','FP','FN','TN','excluded'}, ...
    'Position',[0.855 0.885 0.14 0.03],'Callback',@(o,~)onFilter(o));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','<< Prev', ...
    'Position',[0.855 0.80 0.066 0.05],'Callback',@(o,~)onStep(o,-1));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','Next >>', ...
    'Position',[0.929 0.80 0.066 0.05],'Callback',@(o,~)onStep(o,1));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','Quit', ...
    'Position',[0.855 0.73 0.14 0.05],'ForegroundColor',[0.6 0 0],'Callback',@(o,~)close(ancestor(o,'figure')));
S.hInfo = uicontrol(S.fig,'Style','text','Units','normalized','String','', ...
    'Position',[0.855 0.43 0.14 0.27],'BackgroundColor','w','HorizontalAlignment','left','FontSize',10);

% counts summary
cats = {items.cat};
S.hCounts = uicontrol(S.fig,'Style','text','Units','normalized', ...
    'String',sprintf('TP %d  FP %d\nFN %d  TN %d\nexcl %d', ...
       sum(strcmp(cats,'TP')),sum(strcmp(cats,'FP')),sum(strcmp(cats,'FN')), ...
       sum(strcmp(cats,'TN')),sum(strcmp(cats,'excluded'))), ...
    'Position',[0.855 0.33 0.14 0.09],'BackgroundColor','w','HorizontalAlignment','left','FontSize',10);

% colour legend
legend_str = {
 'LEGEND'
 'title green = correct'
 'title red = incorrect'
 'chan green = GT & detected'
 'chan orange = GT, missed'
 'chan red = detected, not GT'
 'red line = stim on'
 'magenta = stim off'
 'green dash = GT AD'
 'blue band = detection'
 'dashed blue = <3 s run'};
uicontrol(S.fig,'Style','text','Units','normalized','String',legend_str, ...
    'Position',[0.855 0.02 0.14 0.30],'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontSize',8);

set(S.fig,'WindowKeyPressFcn',@onKey);
guidata(S.fig,S);
applyFilter(S.fig, start_filter);
end


function applyFilter(fig, name)
S = guidata(fig);
cats = {S.items.cat};
if strcmp(name,'All'), S.order = 1:numel(S.items);
else, S.order = find(strcmp(cats,name)); end
S.pos = 1;
% sync popup
opts = get(S.hFilter,'String'); k = find(strcmp(opts,name),1);
if ~isempty(k), set(S.hFilter,'Value',k); end
guidata(fig,S);
redraw(fig);
end

function onFilter(o)
fig = ancestor(o,'figure'); opts = get(o,'String');
applyFilter(fig, opts{get(o,'Value')});
end

function onStep(o, d)
fig = ancestor(o,'figure'); S = guidata(fig);
if isempty(S.order), return; end
S.pos = max(1, min(numel(S.order), S.pos + d));
guidata(fig,S); redraw(fig);
end

function onKey(fig, ev)
switch ev.Key
    case 'leftarrow',  onStepFig(fig,-1);
    case 'rightarrow', onStepFig(fig, 1);
    case 'uparrow',    S=guidata(fig); S.gain=S.gain*1.3; guidata(fig,S); redraw(fig);
    case 'downarrow',  S=guidata(fig); S.gain=S.gain/1.3; guidata(fig,S); redraw(fig);
end
end
function onStepFig(fig,d)
S = guidata(fig); if isempty(S.order), return; end
S.pos = max(1,min(numel(S.order),S.pos+d)); guidata(fig,S); redraw(fig);
end


function redraw(fig)
S = guidata(fig);
cla(S.ax);
if isempty(S.order)
    title(S.ax,'(no clips in this category)'); set(S.hInfo,'String',''); return;
end
it = S.items(S.order(S.pos));

clip = load_clip(S, it.clip_name);
fs = clip.fs; labels = clip.labels(:); nCh = numel(labels); V = clip.values;
[bp,bidx] = find_bipolar_pairs(labels,1:nCh);
nP = size(bp,1);

gt_chans = {}; if ~isempty(it.gt_chans), gt_chans = strsplit(it.gt_chans,';'); end

tt = clip.clip_start + (0:size(V,1)-1)/fs;
hold(S.ax,'on');
stim_lab = clip.stim_pair;
% keep only NON-stim bipolar pairs (pairs not sharing a stim contact)
keep = false(nP,1);
for r=1:nP, keep(r) = ~any(ismember(bp(r,:), stim_lab)); end
kp = find(keep);
nShow = numel(kp);
plabs = cell(nShow,1); traces = nan(size(V,1),nShow);
for ii=1:nShow
    r = kp(ii);
    traces(:,ii) = V(:,bidx(r,1)) - V(:,bidx(r,2));
    plabs{ii} = [bp{r,1} '-' bp{r,2}];
end
traces = traces - mean(traces,1,'omitnan');
try, traces = bandstop(traces,[58 62],fs); catch, end
sd = median(std(traces,0,1,'omitnan'),'omitnan'); if sd==0||isnan(sd), sd=1; end
step = 6*sd;
det_here = ~isempty(it.chan);
for c=1:nShow
    in_gt  = ismember(plabs{c}, gt_chans);
    is_det = det_here && strcmp(plabs{c}, it.chan);
    if in_gt && is_det,       col=[0 0.55 0];   % AD channel: GT & detected (agree)
    elseif in_gt && ~is_det,  col=[0.9 0.5 0];  % AD channel: GT only (missed)
    elseif ~in_gt && is_det,  col=[0.85 0 0];   % AD channel: detected only (false)
    else,                     col=[0 0 0]; end
    lw = 0.5 + 1.0*(is_det || in_gt);
    plot(S.ax, tt, traces(:,c)*S.gain - (c-1)*step, 'Color',col,'LineWidth',lw);
    text(S.ax, tt(1), -(c-1)*step, [plabs{c} ' '], 'HorizontalAlignment','right', ...
        'FontSize',9,'Color',col);
end
yl=[-nShow*step, step]; ylim(S.ax,yl); xlim(S.ax,[tt(1) tt(end)]);
% stim markers
plot(S.ax,[clip.stimOn clip.stimOn],yl,'Color',[0.85 0 0],'LineWidth',1.2);
plot(S.ax,[clip.stimOff clip.stimOff],yl,'Color',[0.8 0 0.8],'LineWidth',1.2);
% ground-truth AD interval (green dashed)
if it.gtAD && ~isnan(it.AD_onset)
    plot(S.ax,[it.AD_onset it.AD_onset],yl,'--','Color',[0 0.55 0],'LineWidth',1.5);
    if ~isnan(it.AD_offset)
        plot(S.ax,[it.AD_offset it.AD_offset],yl,'--','Color',[0 0.55 0],'LineWidth',1.0);
    end
    text(S.ax,it.AD_onset,yl(1),' GT AD','Color',[0 0.55 0],'VerticalAlignment','bottom');
end
% near-misses: supra-threshold runs shorter than N windows (<3 s), dashed faint
for q = 1:numel(it.shorts)
    sh = it.shorts(q);
    patch(S.ax,[sh.t0 sh.t1 sh.t1 sh.t0],[yl(1) yl(1) yl(2) yl(2)], ...
        [0.3 0.6 1],'FaceAlpha',0.06,'EdgeColor',[0.3 0.6 1],'LineStyle',':');
    text(S.ax, sh.t0, yl(1), sprintf(' <3s %s', sh.chan), 'Color',[0.2 0.45 0.9], ...
        'FontSize',8,'VerticalAlignment','bottom');
end
% detector: shaded above-threshold window (blue) + onset line
if ~isnan(it.det_t0)
    patch(S.ax,[it.det_t0 it.det_t1 it.det_t1 it.det_t0],[yl(1) yl(1) yl(2) yl(2)], ...
        [0 0 0.85],'FaceAlpha',0.10,'EdgeColor','none');
end
if ~isnan(it.onset)
    plot(S.ax,[it.onset it.onset],yl,'Color',[0 0 0.85],'LineWidth',1.5);
    text(S.ax,it.onset,yl(2),sprintf(' detected %s',it.chan),'Color',[0 0 0.85], ...
        'VerticalAlignment','top');
end
set(S.ax,'YTick',[]); S.ax.XAxis.Exponent=0; xlabel(S.ax,'Time (s)');
[ctxt, ccol] = cat_display(it.cat);
title(S.ax, sprintf('%s   %s   (%d/%d in filter)', it.clip_name, ctxt, S.pos, numel(S.order)), ...
    'Color',ccol,'Interpreter','none');

% info panel
set(S.hInfo,'String',sprintf(['category: %s\nscore: %.2f  (T=%.2f)\n', ...
    'GT AD: %d  ongoing: %d\nGT chans: %s\ndetected: %s\nstim: %s-%s'], ...
    it.cat, it.score, S.p.T, it.gtAD, it.ad_ongoing, ...
    ternary(isempty(it.gt_chans),'-',it.gt_chans), ...
    ternary(isnan(it.onset),'(none)',sprintf('%s  %.1f-%.1f s',it.chan,it.det_t0,it.det_t1)), ...
    stim_lab{1}, stim_lab{2}));
end


function clip = load_clip(S, cname)
if S.remote
    lp = fullfile(S.cache,[cname '.mat']);
    if ~exist(lp,'file')
        [host,rpath] = split_remote(S.remote_clips);
        system(sprintf('scp -q "%s:%s/%s.mat" "%s"', host, rpath, cname, lp));
    end
else
    lp = fullfile(S.clipdir,[cname '.mat']);
end
C = load(lp); clip = C.clip;
end

function [txt, c] = cat_display(cat)
% Title text + color: correct (TP/TN) green, incorrect (FP/FN) red.
switch cat
    case 'TP', txt='TRUE POSITIVE';  c=[0 0.55 0];
    case 'TN', txt='TRUE NEGATIVE';  c=[0 0.55 0];
    case 'FP', txt='FALSE POSITIVE'; c=[0.85 0 0];
    case 'FN', txt='FALSE NEGATIVE'; c=[0.85 0 0];
    otherwise, txt='EXCLUDED';       c=[0.5 0.5 0.5];
end
end

function s = ternary(cond,a,b), if cond, s=a; else, s=b; end, end

function [host,rpath]=split_remote(remote)
k=strfind(remote,':'); k=k(1); host=remote(1:k-1); rpath=regexprep(remote(k+1:end),'/+$','');
end

function rr = true_runs(mask)
% Start/end indices of each maximal run of true values in a logical vector.
mask = mask(:)';
d = diff([false mask false]);
starts = find(d==1);
ends   = find(d==-1) - 1;
rr = [starts(:) ends(:)];
end
