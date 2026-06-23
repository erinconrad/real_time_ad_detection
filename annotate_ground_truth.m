function annotate_ground_truth(overwrite, annotator)
% ANNOTATE_GROUND_TRUTH  Human ground-truth GUI for afterdischarges, one
% CLIP at a time (each clip = one stimulation = one trial).
%
%   annotate_ground_truth                 % annotate, skip clips already done
%   annotate_ground_truth(1)              % re-annotate everything
%   annotate_ground_truth(0,'EC')         % tag annotator initials
%
% For each clip in clips/ it shows the stim-electrode EEG (bipolar,
% notch-filtered) with stim on/off marked, and asks whether an afterdischarge
% (AD) occurred. If yes, you click the AD ONSET then OFFSET on the plot. This
% detector-independent truth is what makes sensitivity/specificity possible.
%
% Output: ground_truth/clip_gt.csv with one row per clip:
%   clip_name, ieeg_name, stim_pair, stimOn, stimOff, AD ('y'/'n'/'u'),
%   AD_onset, AD_offset, Annotator, notes
%
% Controls: a dialog asks Yes / No / Skip-Quit. After "Yes" you click twice
% (onset, offset). "Skip" leaves the clip for later; "Quit" stops and saves.

if nargin < 1 || isempty(overwrite), overwrite = 0; end
if nargin < 2 || isempty(annotator), annotator = getenv('USER'); end

L = rt_paths;

listing = dir(fullfile(L.clip_dir, '*.mat'));
assert(~isempty(listing), 'No clips in %s. Run export_stim_clips first.', L.clip_dir);

gt_csv = fullfile(L.gt_dir, 'clip_gt.csv');
gt = init_or_load_gt(gt_csv);

clip_names = erase({listing.name}, '.mat');

for i = 1:numel(listing)
    cname = clip_names{i};

    % already annotated?
    row = find(strcmp(gt.clip_name, cname), 1);
    if ~isempty(row) && ~overwrite && ismember(gt.AD{row}, {'y','n','u'})
        continue;
    end

    C = load(fullfile(L.clip_dir, listing(i).name)); clip = C.clip;
    fs = clip.fs;
    labels = clip.labels;
    nCh = numel(labels);

    % bipolar montage within the saved electrode contacts
    [~,~,alt] = find_bipolar_pairs(labels, 1:nCh);
    V = clip.values;
    Vb = nan(size(V));
    valid = ~isnan(alt);
    Vb(:,valid) = V(:,valid) - V(:,alt(valid));
    Vb = Vb - mean(Vb,1,'omitnan');
    show = find(valid);
    if isempty(show), Vb = V - mean(V,1,'omitnan'); show = (1:nCh)'; end

    tt = clip.clip_start + (0:size(V,1)-1)/fs;
    Vp = Vb(:, show);
    try, Vp = bandstop(Vp,[58 62],fs); catch, end

    fig = figure('Name',sprintf('%s  (%s-%s)',cname,clip.stim_pair{1},clip.stim_pair{2}), ...
        'Position',[40 40 1500 850]);
    ax = axes(fig); hold(ax,'on');
    sd = nanstd(Vp(:)); if sd==0 || isnan(sd), sd = 1; end
    step = 6*sd;
    for c = 1:numel(show)
        plot(ax, tt, Vp(:,c) - (c-1)*step, 'k');
        text(ax, tt(1), -(c-1)*step, [labels{show(c)} ' '], ...
            'HorizontalAlignment','right','FontSize',10);
    end
    yl = [-(numel(show))*step, step]; ylim(ax, yl);
    plot(ax,[clip.stimOn clip.stimOn],yl,'b','LineWidth',1.5);
    plot(ax,[clip.stimOff clip.stimOff],yl,'r','LineWidth',1.5);
    text(ax,clip.stimOn,yl(2),' stim on','Color','b');
    text(ax,clip.stimOff,yl(2),' stim off','Color','r');
    xlabel(ax,'Time (s)'); set(ax,'YTick',[]); ax.XAxis.Exponent = 0;
    title(ax,sprintf('%s   stim %s-%s   (%d/%d)', cname, ...
        clip.stim_pair{1}, clip.stim_pair{2}, i, numel(listing)),'Interpreter','none');

    choice = questdlg('Afterdischarge present?', sprintf('Clip %d/%d',i,numel(listing)), ...
        'Yes','No','Skip/Quit','No');

    ad = ''; on = NaN; off = NaN; quit_now = false;
    switch choice
        case 'Yes'
            title(ax,'Click AD ONSET, then AD OFFSET');
            [xc,~] = ginput(2); xc = sort(xc);
            ad = 'y'; on = xc(1); off = xc(2);
        case 'No'
            ad = 'n';
        otherwise
            sub = questdlg('Skip this clip or Quit?','Skip/Quit','Skip','Unsure','Quit','Skip');
            if strcmp(sub,'Quit'), quit_now = true;
            elseif strcmp(sub,'Unsure'), ad = 'u';
            end
    end
    if ishandle(fig), close(fig); end

    if ~isempty(ad)
        newrow = {cname, clip.ieeg_name, strjoin(clip.stim_pair,'-'), ...
            clip.stimOn, clip.stimOff, ad, on, off, annotator, ''};
        if isempty(row)
            gt = [gt; newrow]; %#ok<AGROW>
        else
            gt(row,:) = newrow;
        end
        writetable(gt, gt_csv);
    end

    if quit_now
        fprintf('Stopped. Progress saved to %s\n', gt_csv);
        return;
    end
end

fprintf('Done. Ground truth -> %s\n', gt_csv);

end


function gt = init_or_load_gt(gt_csv)
if exist(gt_csv,'file')
    gt = readtable(gt_csv,'TextType','char');
    if ~iscell(gt.AD),        gt.AD        = cellstr(string(gt.AD));        end
    if ~iscell(gt.Annotator), gt.Annotator = cellstr(string(gt.Annotator)); end
    return;
end
gt = table('Size',[0 10], ...
    'VariableTypes',{'cell','cell','cell','double','double','cell','double','double','cell','cell'}, ...
    'VariableNames',{'clip_name','ieeg_name','stim_pair','stimOn','stimOff', ...
                     'AD','AD_onset','AD_offset','Annotator','notes'});
end
