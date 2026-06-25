function annotate_ground_truth(overwrite, annotator, remote_clips)
% ANNOTATE_GROUND_TRUTH  Human ground-truth GUI for afterdischarges, one
% CLIP at a time (each clip = one stimulation = one trial).
%
%   annotate_ground_truth                 % annotate clips in the local clips/ dir
%   annotate_ground_truth(1)              % re-annotate everything
%   annotate_ground_truth(0,'EC')         % tag annotator initials
%   annotate_ground_truth(0,'EC', 'user@host:/path/to/clips')   % REMOTE mode
%
% LOCAL mode (default): reads clip .mat files from the project clips/ folder.
% REMOTE mode: fetches each clip's .mat over scp one at a time, annotates,
% saves to the LOCAL ground_truth/clip_gt.csv, deletes the temp copy. Only the
% .mat is transferred. Use key-based SSH so scp does not prompt for a password.
%
% DISPLAY (all bipolar):
%   * Shows all electrode bipolar pairs by default (uncheck "Show all
%     electrode channels" to hide the stim-artifact pairs). The assigned stim
%     pair is drawn in RED.
%   * Up / Down arrows increase / decrease the gain.
%
% DECISIONS: Yes (then click AD onset, then offset) / No / Unsure / Skip / Quit.
%   * Before clicking Yes, select the channel(s) showing the AD in the
%     "AD channels" list (multi-select: ctrl/shift-click).
%   * "Flag stim-detection problem" checkbox + free-text Notes for later review.
%
% Output: ground_truth/clip_gt.csv, one row per clip:
%   clip_name, ieeg_name, stim_pair, stimOn, stimOff, AD ('y'/'n'/'u'),
%   AD_onset, AD_offset, AD_channels, stim_problem (0/1), Annotator, notes

if nargin < 1 || isempty(overwrite), overwrite = 0; end
if nargin < 2 || isempty(annotator), annotator = getenv('USER'); end
if nargin < 3, remote_clips = ''; end

L = rt_paths;

if isempty(remote_clips)
    try
        cfg = rt_local_config;
        if isfield(cfg,'remote_clips'), remote_clips = cfg.remote_clips; end
    catch
    end
end
remote = ~isempty(remote_clips);

gt_csv = fullfile(L.gt_dir, 'clip_gt.csv');
gt = init_or_load_gt(gt_csv);

if remote
    names = list_remote_clips(remote_clips);
    cache = fullfile(tempdir, 'ad_clip_cache');
    if ~exist(cache,'dir'), mkdir(cache); end
    fprintf('Remote mode: %d clips at %s\n', numel(names), remote_clips);
else
    listing = dir(fullfile(L.clip_dir, '*.mat'));
    assert(~isempty(listing), 'No clips in %s. Run export_stim_clips first.', L.clip_dir);
    names = erase({listing.name}, '.mat');
end

for i = 1:numel(names)
    cname = names{i};

    row = find(strcmp(gt.clip_name, cname), 1);
    if ~isempty(row) && ~overwrite && ismember(gt.AD{row}, {'y','n','u'})
        continue;
    end

    % existing values (so re-shown clips keep flag/notes/channels)
    def = struct('flag',false,'notes','','chans','','ongoing',false);
    if ~isempty(row)
        def.flag    = logical(gt.stim_problem(row));
        def.notes   = gt.notes{row};
        def.chans   = gt.AD_channels{row};
        def.ongoing = logical(gt.ad_ongoing(row));
    end

    if remote
        local_mat = fullfile(cache, [cname '.mat']);
        if ~scp_fetch(remote_clips, [cname '.mat'], local_mat)
            warning('Could not fetch %s, skipping.', cname);
            continue;
        end
    else
        local_mat = fullfile(L.clip_dir, [cname '.mat']);
    end

    C = load(local_mat); clip = C.clip;

    R = annotate_one_clip(clip, cname, i, numel(names), def);

    if ~isempty(R.ad) || R.flag || ~isempty(R.notes)
        newrow = {cname, clip.ieeg_name, strjoin(clip.stim_pair,'-'), ...
            clip.stimOn, clip.stimOff, R.ad, R.on, R.off, double(R.ad_ongoing), ...
            R.ad_channels, double(R.flag), annotator, R.notes};
        if isempty(row), gt = [gt; newrow]; else, gt(row,:) = newrow; end %#ok<AGROW>
        writetable(gt, gt_csv);
    end

    if remote && exist(local_mat,'file'), delete(local_mat); end

    if R.quit
        fprintf('Stopped. Progress saved to %s\n', gt_csv);
        return;
    end
end

fprintf('Done. Ground truth -> %s\n', gt_csv);

end


% ========================================================================
function R = annotate_one_clip(clip, cname, idx, ntot, def)
% Interactive single-clip annotator. Returns struct R with fields:
%   ad, on, off, ad_channels, flag, notes, quit
R = struct('ad','', 'on',NaN, 'off',NaN, 'ad_ongoing',logical(def.ongoing), ...
    'ad_channels','', 'flag',logical(def.flag), 'notes',def.notes, 'quit',false);

fs = clip.fs;
labels = clip.labels(:);
nCh = numel(labels);
V = clip.values;

% --- bipolar pairs on the saved electrode ---
[bp, bidx] = find_bipolar_pairs(labels, 1:nCh);
nP = size(bp,1);
if nP == 0
    traces = V - mean(V,1,'omitnan');
    plabels = labels;
    dmask = true(nCh,1);
    is_stim = ismember(labels, clip.stim_pair);
else
    traces = nan(size(V,1), nP);
    plabels = cell(nP,1);
    is_stim = false(nP,1);
    for r = 1:nP
        traces(:,r) = V(:,bidx(r,1)) - V(:,bidx(r,2));
        plabels{r} = [bp{r,1} '-' bp{r,2}];
        is_stim(r) = all(ismember(bp(r,:), clip.stim_pair));
    end
    dmask = false(nP,1);
    for r = 1:nP
        dmask(r) = ~any(ismember(bp(r,:), clip.stim_pair));
    end
    if ~any(dmask), dmask(:) = true; end
end
traces = traces - mean(traces,1,'omitnan');
try, traces = bandstop(traces,[58 62],fs); catch, end

base_std = median(std(traces(:,dmask),0,1,'omitnan'),'omitnan');
if base_std == 0 || isnan(base_std), base_std = 1; end
tt = clip.clip_start + (0:size(V,1)-1)/fs;

% prefill AD-channel selection
selidx = [];
if ~isempty(def.chans)
    want = strsplit(def.chans, ';');
    selidx = find(ismember(plabels, want));
end

% --- figure + state ---
S.fig = figure('Name',sprintf('%s  (%s-%s)',cname,clip.stim_pair{1},clip.stim_pair{2}), ...
    'Position',[40 40 1500 850], 'Color','w');
S.ax = axes('Parent',S.fig, 'Position',[0.06 0.10 0.78 0.82]); hold(S.ax,'on');
S.tt = tt; S.traces = traces; S.plabels = {plabels}; S.dmask = dmask; S.is_stim = is_stim;
S.step0 = 6*base_std; S.gain = 1; S.show_all = true; S.flag = logical(def.flag);
S.ongoing = logical(def.ongoing);
S.stimOn = clip.stimOn; S.stimOff = clip.stimOff;
S.titlestr = sprintf('%s   stim %s-%s   (%d/%d)', cname, clip.stim_pair{1}, ...
    clip.stim_pair{2}, idx, ntot);
S.result = R;
guidata(S.fig, S);

% --- controls (right panel) ---
uicontrol(S.fig,'Style','checkbox','Units','normalized','String','Show all electrode channels', ...
    'Value',double(S.show_all),'Position',[0.855 0.955 0.14 0.025],'BackgroundColor','w','Callback',@(o,~)onShowAll(o));
uicontrol(S.fig,'Style','checkbox','Units','normalized','String','Flag stim-detection problem', ...
    'Value',double(S.flag),'Position',[0.855 0.927 0.14 0.024],'BackgroundColor','w','Callback',@(o,~)onFlag(o));
uicontrol(S.fig,'Style','checkbox','Units','normalized','String','AD continues past clip end', ...
    'Value',double(S.ongoing),'Position',[0.855 0.899 0.14 0.024],'BackgroundColor','w','Callback',@(o,~)onOngoing(o));
uicontrol(S.fig,'Style','text','Units','normalized','String','Up/Down arrows: gain', ...
    'Position',[0.855 0.872 0.14 0.020],'BackgroundColor','w','HorizontalAlignment','left');
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','Yes - mark AD', ...
    'Position',[0.855 0.815 0.14 0.05],'FontWeight','bold','Callback',@(o,~)onYes(o));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','No AD', ...
    'Position',[0.855 0.762 0.14 0.048],'Callback',@(o,~)onDecide(o,'n'));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','Unsure', ...
    'Position',[0.855 0.710 0.14 0.048],'Callback',@(o,~)onDecide(o,'u'));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','Skip', ...
    'Position',[0.855 0.658 0.14 0.048],'Callback',@(o,~)onDecide(o,''));
uicontrol(S.fig,'Style','pushbutton','Units','normalized','String','Quit', ...
    'Position',[0.855 0.606 0.14 0.048],'ForegroundColor',[0.6 0 0],'Callback',@(o,~)onQuit(o));

uicontrol(S.fig,'Style','text','Units','normalized','String','AD channels (multi-select):', ...
    'Position',[0.855 0.580 0.14 0.020],'BackgroundColor','w','HorizontalAlignment','left');
S.hChans = uicontrol(S.fig,'Style','listbox','Units','normalized','Max',2,'Min',0, ...
    'String',plabels,'Value',selidx,'Position',[0.855 0.345 0.14 0.230]);

uicontrol(S.fig,'Style','text','Units','normalized','String','Notes:', ...
    'Position',[0.855 0.305 0.14 0.022],'BackgroundColor','w','HorizontalAlignment','left');
S.hNotes = uicontrol(S.fig,'Style','edit','Units','normalized','Max',3,'Min',1, ...
    'HorizontalAlignment','left','String',def.notes,'Position',[0.855 0.10 0.14 0.20]);

guidata(S.fig, S);   % store handles
set(S.fig,'WindowKeyPressFcn',@onKey);

replot(S.fig);
uiwait(S.fig);

if ishandle(S.fig)
    S = guidata(S.fig);
    R = S.result; R.flag = S.flag;
    close(S.fig);
end
end


function replot(fig)
S = guidata(fig);
cla(S.ax); hold(S.ax,'on');
if S.show_all, sel = 1:size(S.traces,2); else, sel = find(S.dmask); end
labs = S.plabels{1};
for c = 1:numel(sel)
    r = sel(c);
    if S.is_stim(r), col = [0.85 0 0]; lw = 1.5; else, col = [0 0 0]; lw = 0.5; end
    plot(S.ax, S.tt, S.traces(:,r)*S.gain - (c-1)*S.step0, 'Color', col, 'LineWidth', lw);
    text(S.ax, S.tt(1), -(c-1)*S.step0, [labs{r} ' '], ...
        'HorizontalAlignment','right','FontSize',10,'Color',col);
end
yl = [-(numel(sel))*S.step0, S.step0];
ylim(S.ax, yl); xlim(S.ax, [S.tt(1) S.tt(end)]);
plot(S.ax,[S.stimOn S.stimOn], yl, 'b','LineWidth',1.5);
plot(S.ax,[S.stimOff S.stimOff], yl, 'Color',[0.8 0 0.8],'LineWidth',1.5);
text(S.ax,S.stimOn,yl(2),' stim on','Color','b');
text(S.ax,S.stimOff,yl(2),' stim off','Color',[0.8 0 0.8]);
set(S.ax,'YTick',[]); S.ax.XAxis.Exponent = 0; xlabel(S.ax,'Time (s)');
title(S.ax, S.titlestr, 'Interpreter','none');
end


function onKey(fig, ev)
S = guidata(fig);
switch ev.Key
    case 'uparrow',   S.gain = S.gain*1.3; guidata(fig,S); replot(fig);
    case 'downarrow', S.gain = S.gain/1.3; guidata(fig,S); replot(fig);
end
end

function onShowAll(o)
fig = ancestor(o,'figure'); S = guidata(fig);
S.show_all = logical(get(o,'Value')); guidata(fig,S); replot(fig);
end

function onFlag(o)
fig = ancestor(o,'figure'); S = guidata(fig);
S.flag = logical(get(o,'Value')); guidata(fig,S);
end

function onYes(o)
fig = ancestor(o,'figure'); S = guidata(fig);
if S.ongoing
    % AD runs past the clip end: only need the onset; offset = clip end
    title(S.ax,'Click AD ONSET (AD continues past clip end; click on the AD channel)','Color',[0 0.5 0]);
    figure(fig); [x,y] = ginput(1);
    S.result.ad = 'y'; S.result.on = x; S.result.off = S.tt(end);
    S.result.ad_ongoing = true;
else
    title(S.ax,'Click AD ONSET, then AD OFFSET (click on the AD channel)','Color',[0 0.5 0]);
    figure(fig); [x,y] = ginput(2); [xs,ord] = sort(x); y = y(ord);
    S.result.ad = 'y'; S.result.on = xs(1); S.result.off = xs(2);
    S.result.ad_ongoing = false;
end
% AD channels: explicit list selection wins; otherwise infer from the trace
% nearest to where you clicked on the y-axis
chans = read_channels(S.hChans, S.plabels{1});
if isempty(chans), chans = nearest_channels(S, y); end
S.result.ad_channels = chans;
S.result.notes = read_notes(S.hNotes);
guidata(fig,S); uiresume(fig);
end

function s = nearest_channels(S, yvals)
if S.show_all, sel = 1:size(S.traces,2); else, sel = find(S.dmask); end
labs = S.plabels{1};
base = -(0:numel(sel)-1) * S.step0;     % y baseline of each displayed trace
idx = [];
for k = 1:numel(yvals)
    [~, c] = min(abs(yvals(k) - base));
    idx(end+1) = sel(c); %#ok<AGROW>
end
idx = unique(idx, 'stable');
s = strjoin(labs(idx), ';');
end

function onDecide(o, ad)
fig = ancestor(o,'figure'); S = guidata(fig);
S.result.ad = ad; S.result.on = NaN; S.result.off = NaN; S.result.ad_ongoing = false;
if strcmp(ad,'y'), S.result.ad_channels = read_channels(S.hChans, S.plabels{1}); end
S.result.notes = read_notes(S.hNotes);
guidata(fig,S); uiresume(fig);
end

function onOngoing(o)
fig = ancestor(o,'figure'); S = guidata(fig);
S.ongoing = logical(get(o,'Value')); guidata(fig,S);
end

function onQuit(o)
fig = ancestor(o,'figure'); S = guidata(fig);
S.result.notes = read_notes(S.hNotes);
S.result.quit = true; guidata(fig,S); uiresume(fig);
end

function s = read_notes(h)
v = get(h,'String');
if isempty(v), s = ''; return; end
s = strtrim(strjoin(cellstr(v)', ' '));
end

function s = read_channels(h, labs)
sel = get(h,'Value');
if isempty(sel), s = ''; return; end
s = strjoin(labs(sel), ';');
end


% ========================================================================
function gt = init_or_load_gt(gt_csv)
canon = {'clip_name','ieeg_name','stim_pair','stimOn','stimOff', ...
         'AD','AD_onset','AD_offset','ad_ongoing','AD_channels','stim_problem','Annotator','notes'};
textcols = {'clip_name','ieeg_name','stim_pair','AD','AD_channels','Annotator','notes'};
numcols  = {'stimOn','stimOff','AD_onset','AD_offset','ad_ongoing','stim_problem'};

if exist(gt_csv,'file')
    gt = readtable(gt_csv,'TextType','char');
    % add any missing columns (older CSVs)
    for k = 1:numel(canon)
        if ~ismember(canon{k}, gt.Properties.VariableNames)
            if ismember(canon{k}, numcols)
                gt.(canon{k}) = nan(height(gt),1);
            else
                gt.(canon{k}) = repmat({''}, height(gt),1);
            end
        end
    end
    % coerce types (empty columns can load as numeric NaN, etc.)
    for k = 1:numel(textcols), gt.(textcols{k}) = tocell(gt.(textcols{k})); end
    for k = 1:numel(numcols),  gt.(numcols{k})  = todouble(gt.(numcols{k})); end
    gt.stim_problem(isnan(gt.stim_problem)) = 0;
    gt.ad_ongoing(isnan(gt.ad_ongoing)) = 0;
    gt = gt(:, canon);
    return;
end

gt = table('Size',[0 numel(canon)], ...
    'VariableTypes',{'cell','cell','cell','double','double','cell','double', ...
                     'double','double','cell','double','cell','cell'}, ...
    'VariableNames', canon);
end

function c = tocell(col)
if iscell(col),        c = col(:);
elseif ischar(col),    c = cellstr(col);
elseif isstring(col),  c = cellstr(col);
elseif isnumeric(col), c = repmat({''}, numel(col), 1);  % empty CSV column
else,                  c = cellstr(string(col));
end
end

function d = todouble(col)
if isnumeric(col),  d = col(:);
elseif iscell(col), d = cellfun(@(x) str2double(string(x)), col(:));
else,               d = double(string(col));
end
end


% ========================================================================
function [host, rpath] = split_remote(remote)
k = strfind(remote, ':');
assert(~isempty(k), 'remote_clips must look like user@host:/path/to/clips');
k = k(1);
host  = remote(1:k-1);
rpath = regexprep(remote(k+1:end), '/+$', '');
end

function names = list_remote_clips(remote)
[host, rpath] = split_remote(remote);
cmd = sprintf('ssh %s "ls -1 %s/*.mat 2>/dev/null"', host, rpath);
[st, out] = system(cmd);
if st ~= 0
    error(['Could not list remote clips at %s.\nCheck SSH access (key-based ' ...
           'login recommended).\n%s'], remote, out);
end
lines = strsplit(strtrim(out), newline);
lines = lines(~cellfun(@isempty, lines));
assert(~isempty(lines), 'No .mat clips found at %s', remote);
names = cell(numel(lines),1);
for i = 1:numel(lines)
    [~, nm] = fileparts(strtrim(lines{i}));
    names{i} = nm;
end
end

function ok = scp_fetch(remote, fname, local_path)
[host, rpath] = split_remote(remote);
src = sprintf('%s:%s/%s', host, rpath, fname);
cmd = sprintf('scp -q "%s" "%s"', src, local_path);
st = system(cmd);
ok = (st == 0) && exist(local_path,'file') == 2;
end
