function export_sessions(overwrite, also_write_edf)
% EXPORT_SESSIONS  Pull each HFS stim session from ieeg.org ONCE and save it
% locally, so the rest of the pipeline (annotation, detection, real-time
% simulation) runs offline and reproducibly.
%
%   export_sessions               % skip sessions already exported
%   export_sessions(1)            % re-export everything
%   export_sessions(overwrite, 1) % also write an .edf alongside the .mat
%
% For each row of hfs_sessions.csv this saves:
%   data/<ieeg_name>_<modifier>.mat   containing a struct `session` with
%       .values     [nSamples x nChannels] raw EEG (microvolts)
%       .fs         sampling rate (Hz)
%       .chLabels   {nChannels x 1} cleaned channel labels
%       .start_time absolute ieeg start time (s) of values(1,:)
%       .end_time   absolute ieeg end time (s)
%       .ieeg_name, .modifier, .aT (annotation table)
%
% The .mat is the working format used by every other script. The optional
% EDF is just for portability/sharing.

if nargin < 1 || isempty(overwrite),      overwrite = 0;      end
if nargin < 2 || isempty(also_write_edf), also_write_edf = 0; end

L = rt_paths;   % sets up paths (project + IEEG toolbox)
pwfile     = L.ieeg_pw_file;
login_name = L.ieeg_login;

fT = readtable(L.session_csv);
nfiles = size(fT,1);

for i = 1:nfiles
    ieeg_name  = fT.ieeg_name{i};
    modifier   = fT.Modifier(i);
    start_time = fT.start(i);
    % the end column may be named 'end', 'xEnd', 'End', etc. depending on the
    % MATLAB version / naming rule. Look it up by name (dynamic field access
    % works even when the name is the reserved word 'end').
    vn = fT.Properties.VariableNames;
    endcol = vn{find(ismember(lower(vn), {'end','xend'}), 1)};
    end_time = fT.(endcol)(i);

    out_mat = fullfile(L.data_dir, sprintf('%s_%d.mat', ieeg_name, modifier));
    if overwrite == 0 && exist(out_mat,'file')
        fprintf('Already exported %s_%d, skipping.\n', ieeg_name, modifier);
        continue;
    end

    fprintf('Exporting %s_%d (%d of %d): %.1f-%.1f s ...\n', ...
        ieeg_name, modifier, i, nfiles, start_time, end_time);

    % Pull (chunked to respect ieeg.org request limits)
    session = pull_session(ieeg_name, start_time, end_time, modifier);

    save(out_mat, 'session', '-v7.3');
    fprintf('  saved %s\n', out_mat);

    if also_write_edf
        try
            write_session_edf(session, fullfile(L.data_dir, ...
                sprintf('%s_%d.edf', ieeg_name, modifier)));
        catch ME
            fprintf('  [EDF write skipped: %s]\n', ME.message);
        end
    end
end

fprintf('\nDone exporting %d sessions to %s\n', nfiles, L.data_dir);

end


function write_session_edf(session, edf_path)
% Write a session to EDF (requires R2021a+). Channel labels are sanitized.
% NaNs are zero-filled (EDF cannot store NaN). Mainly for portability.
fs = session.fs;
v  = session.values;
v(isnan(v)) = 0;

% sanitize labels to valid EDF names
labs = matlab.lang.makeValidName(session.chLabels);
tt = array2timetable(v, 'SampleRate', fs, 'VariableNames', labs);
edfw = edfwrite(edf_path, edfheader("EDF"), tt); %#ok<NASGU>
fprintf('  wrote %s\n', edf_path);
end
