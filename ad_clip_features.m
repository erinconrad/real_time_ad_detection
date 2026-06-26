function [F, dbg] = ad_clip_features(clip, params)
% AD_CLIP_FEATURES  Compute per-window line-length z-scores for the candidate
% AD channels of one clip. This is the front end shared by the offline
% detector and the validation harness.
%
%   F = ad_clip_features(clip, params)
%
% Steps (all causal, so identical online and offline):
%   1. pick the nearest non-stim bipolar pairs (ad_candidate_pairs)
%   2. causal 1-50 Hz band-pass
%   3. robust pre-stim baseline per channel (median/MAD over trimmed windows)
%   4. line length in win_s windows post-stim, z-scored to baseline
%   5. high-frequency artifact guard marks contaminated windows invalid
%
% OUTPUT F:
%   .chan_labels {1xm}      candidate bipolar pair labels
%   .win_t       [1xnWin]   absolute start time of each post-stim window (s)
%   .z           [nWin x m] LL z-score vs pre-stim baseline
%   .valid       [nWin x m] logical, false where HF artifact vetoes the window

p = ad_params(params);
fs = clip.fs;
labels = clip.labels(:);
V = clip.values;
dbg = struct();

[pairs, plabs] = ad_candidate_pairs(labels, clip.stim_pair, p.n_candidates);
F.chan_labels = plabs(:)';
m = size(pairs,1);
if m == 0
    F.win_t = []; F.z = []; F.valid = []; return;
end

% bipolar candidate signals
X = zeros(size(V,1), m);
for j = 1:m, X(:,j) = V(:,pairs(j,1)) - V(:,pairs(j,2)); end
X(isnan(X)) = 0;
Xraw = X;

% causal band-pass for line length
[bb, aa] = butter(p.bp_order, p.bp_band/(fs/2), 'bandpass');
Xf = filter(bb, aa, X);

% causal 60 Hz notch so line noise does not inflate the line length
if p.notch60 && (fs/2) > 62
    [bn, an] = butter(2, [58 62]/(fs/2), 'stop');
    Xf = filter(bn, an, Xf);
end

% sample indices for stim on/off relative to clip start
nS = size(V,1);
sonIdx  = max(1, round((clip.stimOn  - clip.clip_start)*fs) + 1);
soffIdx = max(1, round((clip.stimOff - clip.clip_start)*fs) + 1);

winN = round(p.win_s*fs);

% window start indices: baseline (pre stim onset), post (from stim offset)
bl_starts   = 1:winN:(sonIdx - winN);
post_starts = soffIdx:winN:(nS - winN + 1);
if isempty(post_starts), F.win_t = []; F.z = []; F.valid = []; return; end

% HF guard band (skip if fs too low)
hf_hi = min(p.hf_band(2), floor(fs/2) - 5);
do_hf = p.hf_guard && hf_hi > p.hf_band(1) + 5;

nWin = numel(post_starts);
F.z = nan(nWin, m);
F.valid = true(nWin, m);
F.win_t = clip.clip_start + (post_starts - 1)/fs;
hf_bad   = false(nWin, m);   % windows vetoed by the HF guard
line_bad = false(1, m);      % channels excluded for 60 Hz line noise

% --- next-stim cap (preferred): truncate analysis before the NEXT stim ---
% next_stim_on is the onset time of the following stim in the session, taken
% from the stim event list (set by build_ad_features / export_stim_clips).
post_ok = true(nWin,1);
nso = inf;
if isfield(clip,'next_stim_on') && ~isempty(clip.next_stim_on) && isfinite(clip.next_stim_on)
    nso = clip.next_stim_on;
end
if isfinite(nso)
    post_ok(F.win_t >= (nso - p.guard_margin_s)) = false;
end

% --- optional amplitude-based guard (off unless stim_sat_frac is finite) ---
if isfinite(p.stim_sat_frac)
    sidx = find(ismember(labels, clip.stim_pair));
    if ~isempty(sidx)
        A = max(max(abs(V(sonIdx:min(soffIdx,nS), sidx)), [], 1), [], 2);
        if ~isempty(A) && A > 0
            over = false(nWin,1);
            for w = 1:nWin
                seg = V(post_starts(w):post_starts(w)+winN-1, sidx);
                over(w) = max(abs(seg(:))) > p.stim_sat_frac * A;
            end
            deadwin = ceil(p.guard_deadtime_s / p.win_s);
            if deadwin >= 1, over(1:min(deadwin,nWin)) = false; end
            cutoff = [];
            for w = 1:nWin-1
                if over(w) && over(w+1), cutoff = w; break; end
            end
            if ~isempty(cutoff), post_ok(cutoff:end) = false; end
        end
    end
end

for j = 1:m
    % --- baseline LL (robust, trimmed) ---
    blLL = window_LL(Xf(:,j), bl_starts, winN);
    if isempty(blLL)
        med = 0; scale = 1;
    else
        s = sort(blLL);
        keepN = max(1, floor(numel(s)*(1 - p.baseline_trim)));
        s = s(1:keepN);
        med = median(s);
        scale = 1.4826*median(abs(s - med)) + eps;
    end
    postLL = window_LL(Xf(:,j), post_starts, winN);
    F.z(:,j) = (postLL - med)/scale;

    % --- HF buzz guard: veto windows where HF power dominates the in-band
    % power (buzz), but NOT spiky ADs whose energy is mostly 1-50 Hz ---
    if do_hf
        hfp    = window_bp(Xraw(:,j), post_starts, winN, fs, [p.hf_band(1) hf_hi]);
        inband = window_bp(Xraw(:,j), post_starts, winN, fs, [1 min(50,floor(fs/2)-1)]);
        hf_bad(:,j) = hfp ./ (inband + eps) > p.hf_ratio;
    end
    % exclude this candidate channel entirely if it is swamped by 60 Hz line noise
    if isfinite(p.line_ratio) && (fs/2) > 61
        p60  = bandpower(Xraw(:,j), fs, [59 61]);
        psig = bandpower(Xraw(:,j), fs, [1 min(50, floor(fs/2)-1)]);
        line_bad(j) = p60/(psig + eps) > p.line_ratio;
    end

    F.valid(:,j) = ~hf_bad(:,j) & ~line_bad(j) & post_ok;
end

% diagnostics (why windows were invalid)
dbg.chan_labels = F.chan_labels;
dbg.win_t = F.win_t;
dbg.z = F.z;
dbg.hf_bad = hf_bad;
dbg.line_bad = line_bad;
dbg.cap_bad = ~post_ok;
dbg.T = p.T; dbg.N = p.N;
end


function v = window_LL(sig, starts, winN)
v = nan(numel(starts),1);
for i = 1:numel(starts)
    seg = sig(starts(i):starts(i)+winN-1);
    v(i) = sum(abs(diff(seg)));
end
end

function v = window_bp(sig, starts, winN, fs, band)
v = nan(numel(starts),1);
for i = 1:numel(starts)
    seg = sig(starts(i):starts(i)+winN-1);
    v(i) = bandpower(seg, fs, band);
end
end
