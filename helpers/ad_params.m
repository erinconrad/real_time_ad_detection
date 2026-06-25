function p = ad_params(o)
% AD_PARAMS  Default parameters for the line-length AD detector, overridable
% by passing a struct o with any subset of fields.

d.bp_band       = [1 50];   % causal band-pass before line length (Hz)
d.bp_order      = 4;        % Butterworth order
d.win_s         = 1;        % LL window length (s)
d.n_candidates  = 2;        % nearest non-stim bipolar pairs to search
d.baseline_trim = 0.2;      % drop this top fraction of pre-stim windows (artifact)
d.hf_guard      = true;     % veto windows with high high-frequency power
d.hf_band       = [70 120]; % HF artifact band (Hz)
d.hf_z          = 5;        % HF z above baseline => window invalid
d.N             = 3;        % consecutive windows above threshold to call AD
d.T             = 4;        % z-score threshold (tuned by validate_ad_detector)
d.min_ad_dur    = 3;        % annotated ADs shorter than this excluded from validation (s)

if nargin < 1 || isempty(o), o = struct(); end
p = d;
f = fieldnames(o);
for i = 1:numel(f), p.(f{i}) = o.(f{i}); end
end
