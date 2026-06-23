function clipTable = test_export_clips(session_row, n_clips, overwrite, chmode, params)
% TEST_EXPORT_CLIPS  Export just a few clips from one session and plot each,
% so you can confirm saving + plotting before running the full export.
%
%   test_export_clips                 % session 1, first 2 clips
%   test_export_clips(1, 3)           % session 1, first 3 clips
%   test_export_clips(row, n, overwrite, chmode, params)
%
% It calls export_stim_clips with a max_clips cap, then opens plot_clip on each
% clip it just saved. Re-run the full export later with export_stim_clips.

if nargin < 1 || isempty(session_row), session_row = 1; end
if nargin < 2 || isempty(n_clips),     n_clips = 2;     end
if nargin < 3 || isempty(overwrite),   overwrite = 1;   end  % default: re-make during testing
if nargin < 4 || isempty(chmode),      chmode = 'electrode'; end
if nargin < 5, params = struct(); end

clipTable = export_stim_clips(session_row, [], [], overwrite, chmode, params, n_clips);

if isempty(clipTable)
    warning('No clips were exported (no stims found?). Try test_stim_detection first.');
    return;
end

fprintf('\nPlotting %d exported clip(s)...\n', height(clipTable));
for i = 1:height(clipTable)
    plot_clip(clipTable.clip_name{i});
end

end
