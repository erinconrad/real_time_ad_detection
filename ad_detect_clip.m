function out = ad_detect_clip(clip, params)
% AD_DETECT_CLIP  Line-length afterdischarge detector for one clip.
%
%   out = ad_detect_clip(clip, params)
%
% Declares an AD if any candidate channel has N consecutive post-stim windows
% whose line-length z-score (vs pre-stim baseline) exceeds threshold T.
% Defaults (ad_params): 1-50 Hz band-pass, 1-s windows, N=3, T=4, nearest 2
% non-stim bipolar pairs, HF artifact guard.
%
% OUTPUT out:
%   .ad       logical  (AD detected)
%   .onset    double   absolute onset time (s) of the detected run, else NaN
%   .channel  char     channel of the detected AD, else ''
%   .score    double   best run min-z (compare to threshold T)
%   .F        struct   features (for plotting/debug)

p = ad_params(params);
F = ad_clip_features(clip, p);
[score, onset_t, chan] = ad_apply_rule(F, p.N);

out.ad      = score > p.T;
out.score   = score;
out.onset   = NaN; out.channel = '';
if out.ad
    out.onset = onset_t; out.channel = chan;
end
out.F = F;
end
