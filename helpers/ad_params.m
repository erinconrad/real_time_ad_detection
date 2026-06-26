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
d.notch60       = true;     % causal 60 Hz notch before line length
d.line_ratio    = 0.5;      % exclude a candidate channel if its 60 Hz band power
                            % exceeds this fraction of its 1-50 Hz power (Inf=off)
d.guard_margin_s = 0.5;     % cut post-stim analysis this many s BEFORE the next
                            % stim onset (next-stim time comes from the stim
                            % event list, not a signal check)
d.stim_sat_frac = Inf;      % optional amplitude-based next-stim guard; OFF by
                            % default (Inf). Set e.g. 0.8 to also enable it.
d.guard_deadtime_s = 2;     % amplitude guard: ignore this much post-offset decay
d.N             = 3;        % consecutive windows above threshold to call AD
d.T             = 4;        % z-score threshold (tuned by validate_ad_detector)
d.min_ad_dur    = 3;        % annotated ADs shorter than this excluded from validation (s)

if nargin < 1 || isempty(o), o = struct(); end
p = d;
f = fieldnames(o);
for i = 1:numel(f), p.(f{i}) = o.(f{i}); end
end
