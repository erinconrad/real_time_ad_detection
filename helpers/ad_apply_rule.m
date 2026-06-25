function [score, onset_t, chan] = ad_apply_rule(F, N)
% AD_APPLY_RULE  Apply the "N consecutive windows" rule to clip features.
%
%   [score, onset_t, chan] = ad_apply_rule(F, N)
%
% F is the struct from ad_clip_features (fields z [nWin x nChan], valid, win_t,
% chan_labels). For each candidate channel it finds the run of N consecutive
% VALID windows with the highest minimum z-score; `score` is that best
% min-z over all channels (so a clip is an AD at threshold T iff score > T,
% which makes threshold tuning trivial). onset_t is the start time of that
% best run and chan its channel.

score = -inf; onset_t = NaN; chan = '';
[nWin, m] = size(F.z);
if nWin < N, return; end

for c = 1:m
    zc = F.z(:,c);
    zc(~F.valid(:,c)) = -inf;     % invalid (artifact) windows break runs
    for s = 1:(nWin - N + 1)
        rm = min(zc(s:s+N-1));
        if rm > score
            score = rm; onset_t = F.win_t(s); chan = F.chan_labels{c};
        end
    end
end
end
