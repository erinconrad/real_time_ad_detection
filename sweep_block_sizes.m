function T = sweep_block_sizes(block_ms_list)
% SWEEP_BLOCK_SIZES  Run the real-time sim + evaluation across several block
% sizes to show how detection accuracy and latency trade off with how often
% data is delivered (emulating different real-time polling intervals).
%
%   T = sweep_block_sizes              % default [100 250 500 1000] ms
%   T = sweep_block_sizes([100 500])
%
% Returns (and prints) a table with sensitivity, specificity, FPR, onset
% error, and worst-case per-block CPU time for each block size. Good single
% figure/table for the prelim data.

if nargin < 1 || isempty(block_ms_list), block_ms_list = [100 250 500 1000]; end

rows = [];
for k = 1:numel(block_ms_list)
    bms = block_ms_list(k);
    res = run_realtime_sim(bms);
    s   = evaluate_detector(bms);
    rows = [rows; table(bms, 100*s.sensitivity, 100*s.specificity, 100*s.fpr, ...
        s.fp_per_min, s.median_onset_err_s, max([res.max_block_cpu_ms]), ...
        median([res.real_time_factor]), ...
        'VariableNames', {'block_ms','sens_pct','spec_pct','fpr_pct', ...
        'fp_per_min','median_onset_err_s','max_block_cpu_ms','median_rt_factor'})]; %#ok<AGROW>
end
T = rows;

fprintf('\n================ Block-size sweep ================\n');
disp(T);

L = rt_paths;
writetable(T, fullfile(L.results_dir,'block_size_sweep.csv'));

end
