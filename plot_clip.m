function plot_clip(clip)
% PLOT_CLIP  Visualize a saved stim clip to confirm the export looks right.
%
%   plot_clip('HUP212_phaseII_1_ev01_LB3-LB4')   % by clip name (in clips/)
%   plot_clip('/full/path/to/clip.mat')          % by .mat path
%   plot_clip(clipStruct)                         % a loaded `clip` struct
%
% Shows the stim-electrode contacts (bipolar, notch-filtered) stacked, with
% the detected stim on/off marked. Channels/fs/times come from the clip's own
% metadata, so this also confirms the saved file is self-describing.

% --- resolve input into a clip struct ---
if ischar(clip) || isstring(clip)
    clip = char(clip);
    if exist(clip,'file') == 2
        S = load(clip);
    else
        L = rt_paths;
        f = fullfile(L.clip_dir, [erase(clip,'.mat') '.mat']);
        assert(exist(f,'file')==2, 'Clip not found: %s', f);
        S = load(f);
    end
    clip = S.clip;
end

fs = clip.fs;
labels = clip.labels;
nCh = numel(labels);

% bipolar within saved contacts
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

figure('Name',sprintf('%s_%d ev%d  %s-%s', clip.ieeg_name, clip.modifier, ...
    clip.event, clip.stim_pair{1}, clip.stim_pair{2}), 'Position',[60 60 1400 800]);
hold on;
sd = nanstd(Vp(:)); if sd==0 || isnan(sd), sd = 1; end
step = 6*sd;
for c = 1:numel(show)
    plot(tt, Vp(:,c) - (c-1)*step, 'k');
    text(tt(1), -(c-1)*step, [labels{show(c)} ' '], ...
        'HorizontalAlignment','right','FontSize',10);
end
yl = [-(numel(show))*step, step]; ylim(yl);
plot([clip.stimOn clip.stimOn], yl, 'b', 'LineWidth',1.5);
plot([clip.stimOff clip.stimOff], yl, 'r', 'LineWidth',1.5);
text(clip.stimOn, yl(2),' stim on','Color','b');
text(clip.stimOff, yl(2),' stim off','Color','r');
set(gca,'YTick',[]); ax = gca; ax.XAxis.Exponent = 0;
xlabel('Time (s)');
title(sprintf('%s_%d  event %d   stim %s-%s   (%d ch, %.1f s)', ...
    clip.ieeg_name, clip.modifier, clip.event, clip.stim_pair{1}, clip.stim_pair{2}, ...
    numel(labels), clip.clip_end-clip.clip_start), 'Interpreter','none');

end
