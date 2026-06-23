function write_edf_clip(edf_path, clip, fs, labels)
% WRITE_EDF_CLIP  Write a [samples x channels] clip to EDF+ (R2021a+).
%
%   write_edf_clip(edf_path, clip, fs, labels)
%
% EDF stores an integer number of fixed-duration data records and cannot
% store NaN. We therefore: (1) replace NaN with 0, and (2) zero-pad the end
% so the clip length is an integer number of 1-second records. The true
% (unpadded) length should be recorded separately in the clip metadata.

if ~exist('edfwrite','file')
    error('edfwrite not found (needs MATLAB R2021a+ Signal Processing Toolbox).');
end

clip(isnan(clip)) = 0;

% pad to whole seconds
nPad = mod(-size(clip,1), round(fs));
if nPad > 0
    clip = [clip; zeros(nPad, size(clip,2))];
end

% valid, unique EDF signal names
labs = matlab.lang.makeValidName(labels(:)');
labs = matlab.lang.makeUniqueStrings(labs);

tt = array2timetable(clip, 'SampleRate', fs, 'VariableNames', labs);

hdr = edfheader("EDF");
hdr.Patient = 'anon';
hdr.Recording = 'stim_clip';

edfw = edfwrite(edf_path, hdr, tt, 'DataRecordDuration', 1); %#ok<NASGU>

end
