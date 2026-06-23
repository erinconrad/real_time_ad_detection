function [summary, perClip] = evaluate_detector(block_ms, post_window_s)
% EVALUATE_DETECTOR  Compare the real-time detector against human ground
% truth and compute sensitivity, specificity, FPR, PPV, and onset-timing
% error. The trial unit is the clip (one stimulation).
%
%   [summary, perClip] = evaluate_detector            % default 250 ms blocks
%   [summary, perClip] = evaluate_detector(block_ms)
%   [summary, perClip] = evaluate_detector(block_ms, post_window_s)
%
% Requires:
%   ground_truth/clip_gt.csv               (from annotate_ground_truth)
%   results/<clip>_rt_<block_ms>ms.mat     (from run_realtime_sim)
%
% A clip counts as a detector POSITIVE if the detector produced any AD
% detection between stim onset and (stim offset + post_window_s). Ground
% truth positive = AD=='y'. Clips marked 'u' or unannotated are excluded.

if nargin < 1 || isempty(block_ms),      block_ms = 250; end
if nargin < 2 || isempty(post_window_s), post_window_s = 5; end

L = rt_paths;

gt_csv = fullfile(L.gt_dir, 'clip_gt.csv');
assert(exist(gt_csv,'file')==2, 'No ground truth at %s. Run annotate_ground_truth first.', gt_csv);
gt = readtable(gt_csv, 'TextType','char');
if ~iscell(gt.AD), gt.AD = cellstr(string(gt.AD)); end

perClip = table();
total_monitored_min = 0;

for i = 1:height(gt)
    cname = gt.clip_name{i};
    rt_file = fullfile(L.results_dir, sprintf('%s_rt_%dms.mat', cname, block_ms));
    if ~exist(rt_file,'file')
        continue;   % clip not yet run through the detector at this block size
    end
    R = load(rt_file); events = R.events;
    ad = events(strcmp(events.Type,'AD'), :);
    ad_times = ad.OnTime;

    stimOn  = gt.stimOn(i);
    stimOff = gt.stimOff(i);
    win = [stimOn, stimOff + post_window_s];
    total_monitored_min = total_monitored_min + (win(2)-win(1))/60;

    in_win = ad_times >= win(1) & ad_times <= win(2);
    detected = any(in_win);
    det_onset = NaN; if detected, det_onset = min(ad_times(in_win)); end

    perClip = [perClip; table({cname}, gt.AD(i), detected, ...
        gt.AD_onset(i), det_onset, ...
        'VariableNames', {'clip','gt_AD','detected','gt_onset','det_onset'})]; %#ok<AGROW>
end

assert(~isempty(perClip), ...
    'No matching detector results for block=%d ms. Run run_realtime_sim(%d) first.', block_ms, block_ms);

% restrict to scored clips (y/n)
scored = ismember(perClip.gt_AD, {'y','n'});
P = perClip(scored, :);
is_pos = strcmp(P.gt_AD,'y');
det    = P.detected;

TP = sum(is_pos & det);
FN = sum(is_pos & ~det);
FP = sum(~is_pos & det);
TN = sum(~is_pos & ~det);

sens = TP/(TP+FN);
spec = TN/(TN+FP);
ppv  = TP/(TP+FP);
fpr  = FP/(FP+TN);
acc  = (TP+TN)/(TP+TN+FP+FN);
fp_per_min = FP / max(total_monitored_min, eps);

tp_rows = is_pos & det;
onset_err = P.det_onset(tp_rows) - P.gt_onset(tp_rows);
onset_err = onset_err(~isnan(onset_err));

summary = struct('block_ms',block_ms,'n_clips',height(P), ...
    'TP',TP,'FN',FN,'FP',FP,'TN',TN, ...
    'sensitivity',sens,'specificity',spec,'ppv',ppv,'fpr',fpr,'accuracy',acc, ...
    'fp_per_min',fp_per_min,'median_onset_err_s',median(onset_err), ...
    'iqr_onset_err_s',[prctile(onset_err,25) prctile(onset_err,75)]);

fprintf('\n================ Detector evaluation (block = %d ms) ================\n', block_ms);
fprintf('Scored clips: %d   (TP %d, FN %d, FP %d, TN %d)\n', height(P), TP, FN, FP, TN);
fprintf('Sensitivity : %.1f%%\n', 100*sens);
fprintf('Specificity : %.1f%%\n', 100*spec);
fprintf('PPV         : %.1f%%\n', 100*ppv);
fprintf('FPR         : %.1f%%  (%.2f false positives / min monitored)\n', 100*fpr, fp_per_min);
fprintf('Accuracy    : %.1f%%\n', 100*acc);
if ~isempty(onset_err)
    fprintf('AD onset timing error (detector - human): median %.2f s (IQR %.2f to %.2f)\n', ...
        median(onset_err), prctile(onset_err,25), prctile(onset_err,75));
end

save(fullfile(L.results_dir, sprintf('eval_%dms.mat', block_ms)), 'summary','perClip');
writetable(perClip, fullfile(L.results_dir, sprintf('eval_%dms_perClip.csv', block_ms)));

end
