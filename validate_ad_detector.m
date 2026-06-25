function res = validate_ad_detector(params)
% VALIDATE_AD_DETECTOR  Patient-grouped cross-validation of the LL AD detector
% against the human ground truth, using cached features (build_ad_features).
%
%   res = validate_ad_detector            % uses ad_params defaults (N=3)
%   res = validate_ad_detector(params)    % override (e.g. struct('N',3,'min_ad_dur',3))
%
% Trial unit = clip. Exclusions (per your annotation rules):
%   * carryover ADs (AD_onset < stimOn) -- AD belongs to the PRIOR stim
%   * annotated ADs shorter than min_ad_dur (non-ongoing) -- not scored
% Threshold T is tuned by LEAVE-ONE-PATIENT-OUT CV (max Youden's J on the
% training patients), so reported metrics are out-of-sample and not leaked
% across the correlated clips of one electrode.
%
% Reports sensitivity, specificity, PPV, FPR (incl. per-minute), accuracy,
% ROC-AUC, and AD onset-timing error; saves results/ad_validation.mat + a
% per-clip CSV.

if nargin < 1, params = struct(); end
L = rt_paths;            % sets up the path (incl. helpers/) first
p = ad_params(params);

f = fullfile(L.results_dir,'ad_features.mat');
assert(exist(f,'file')==2, 'Run build_ad_features first (no %s).', f);
Sload = load(f); feats = Sload.feats;
n = numel(feats);

% --- per-clip score / onset for the fixed N ---
score = -inf(n,1); onset = nan(n,1);
for i = 1:n
    [sc, on] = ad_apply_rule(feats(i).F, p.N);
    score(i) = sc; onset(i) = on;
end

label   = logical([feats.label]');
patient = string({feats.patient}');

% --- exclusions ---
excl = false(n,1);
for i = 1:n
    if label(i)
        if feats(i).AD_onset < feats(i).stimOn, excl(i) = true; end          % carryover
        dur = feats(i).AD_offset - feats(i).AD_onset;
        if ~feats(i).ad_ongoing && dur < p.min_ad_dur, excl(i) = true; end   % too short
    end
end
eval = ~excl;
fprintf('Evaluable clips: %d (%d AD, %d no-AD); excluded %d.\n', ...
    sum(eval), sum(eval & label), sum(eval & ~label), sum(excl));

% --- leave-one-patient-out CV ---
groups = unique(patient(eval));
pred = false(n,1); used = false(n,1); Tfold = [];
for g = 1:numel(groups)
    te = eval & patient==groups(g);
    tr = eval & patient~=groups(g);
    if ~any(tr & label) || ~any(tr & ~label)
        Tstar = p.T;   % can't tune without both classes; fall back to default
    else
        Tstar = pick_threshold(score(tr), label(tr));
    end
    pred(te) = score(te) > Tstar;
    used(te) = true;
    Tfold(end+1) = Tstar; %#ok<AGROW>
end

% --- metrics over held-out predictions ---
y = label(used); yh = pred(used);
TP = sum(y & yh); FN = sum(y & ~yh); FP = sum(~y & yh); TN = sum(~y & ~yh);
sens = TP/max(TP+FN,1); spec = TN/max(TN+FP,1);
ppv  = TP/max(TP+FP,1); fpr = FP/max(FP+TN,1);
acc  = (TP+TN)/max(TP+TN+FP+FN,1);

% false positives per minute monitored (over evaluable no-AD clips)
mon_min = 0;
ev_idx = find(eval);
for i = ev_idx'
    if ~label(i), mon_min = mon_min + (feats(i).clip_end - feats(i).stimOff)/60; end
end
fp_per_min = FP / max(mon_min, eps);

% onset error on true positives (full-length masks)
te_tp = used & label & pred;
oe = onset(te_tp) - arrayfun(@(s) s.AD_onset, feats(te_tp)');
oe = oe(~isnan(oe));

auc = roc_auc(score(eval), label(eval));

% --- report ---
fprintf('\n============ AD detector validation (LOPO CV, N=%d) ============\n', p.N);
fprintf('TP %d  FN %d  FP %d  TN %d\n', TP, FN, FP, TN);
fprintf('Sensitivity : %.1f%%\n', 100*sens);
fprintf('Specificity : %.1f%%\n', 100*spec);
fprintf('PPV         : %.1f%%\n', 100*ppv);
fprintf('FPR         : %.1f%%  (%.2f false pos / min)\n', 100*fpr, fp_per_min);
fprintf('Accuracy    : %.1f%%\n', 100*acc);
fprintf('ROC-AUC     : %.3f\n', auc);
fprintf('Suggested threshold T (median across folds): %.2f\n', median(Tfold));
if ~isempty(oe)
    fprintf('AD onset error (detector - human): median %.2f s (IQR %.2f to %.2f)\n', ...
        median(oe), prctile(oe,25), prctile(oe,75));
end

res = struct('N',p.N,'TP',TP,'FN',FN,'FP',FP,'TN',TN,'sensitivity',sens, ...
    'specificity',spec,'ppv',ppv,'fpr',fpr,'accuracy',acc,'auc',auc, ...
    'fp_per_min',fp_per_min,'T_suggested',median(Tfold), ...
    'onset_err_median',median(oe));

% per-clip table
T = table({feats.clip_name}', patient, label, eval, score, pred, onset, ...
    'VariableNames',{'clip','patient','gt_AD','evaluable','score','pred','onset'});
writetable(T, fullfile(L.results_dir,'ad_validation_perClip.csv'));
save(fullfile(L.results_dir,'ad_validation.mat'),'res','T','Tfold');
end


function Tstar = pick_threshold(score, label)
% threshold maximizing Youden's J (sens + spec - 1)
cand = unique(score(isfinite(score)));
cand = [min(cand)-1; cand(:)];
P = sum(label); Nn = sum(~label);
bestJ = -inf; Tstar = cand(1);
for t = 1:numel(cand)
    yh = score > cand(t);
    sens = sum(label & yh)/max(P,1);
    spec = sum(~label & ~yh)/max(Nn,1);
    J = sens + spec - 1;
    if J > bestJ, bestJ = J; Tstar = cand(t); end
end
end


function a = roc_auc(score, label)
% AUC via Mann-Whitney U on ranks (handles ties)
P = sum(label); Nn = sum(~label);
if P==0 || Nn==0, a = NaN; return; end
s = score; s(~isfinite(s)) = min(s(isfinite(s)))-1;  % push -inf to bottom
r = tiedrank(s);
a = (sum(r(label)) - P*(P+1)/2) / (P*Nn);
end
