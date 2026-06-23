function run_pipeline(block_ms)
% RUN_PIPELINE  Convenience driver for the real-time AD prelim analysis.
%
% Typical order of operations:
%   1) export_sessions          % once: pull data from ieeg.org -> local .mat
%   2) annotate_ground_truth    % human: mark AD y/n + onset/offset per stim
%   3) run_realtime_sim(block)  % machine: streaming detection + timing
%   4) evaluate_detector(block) % stats: sens/spec/FPR/PPV/timing
%
% run_pipeline(block_ms) runs steps 3-4 (assumes 1-2 already done). It does
% NOT run export or annotation, since those need ieeg.org access and a human.

if nargin < 1 || isempty(block_ms), block_ms = 250; end

run_realtime_sim(block_ms);
evaluate_detector(block_ms);

end
