function write_edf_clip(edf_path, clip, fs, labels)
% WRITE_EDF_CLIP  Write a [samples x channels] clip to an EDF file (R2021a+).
%
%   write_edf_clip(edf_path, clip, fs, labels)
%
% Follows the documented edfwrite API: edfwrite(filename, hdr, sigdata) where
% hdr is a header struct from edfheader and sigdata is a NUMERIC MATRIX
% (samples x channels). We write the whole clip as a single data record whose
% duration is the clip length, which avoids any integer-samples-per-record
% constraint and works for any sampling rate. Signal data are written with
% physical scaling (microvolts) into the int16 digital range.
%
% EDF cannot store NaN, so NaNs are replaced with 0. The lossless copy (with
% NaNs and true fs) lives in the clip .mat.

if ~exist('edfwrite','file')
    error('edfwrite not found (needs MATLAB R2021a+ Signal Processing Toolbox).');
end

data = clip;
data(isnan(data)) = 0;
[nSamp, nCh] = size(data);

% valid, unique EDF signal names (<=16 chars), as a string vector
labs = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(labels(:)'));
labs = string(labs);

% physical range per channel. EDF stores PhysicalMin/Max as 8-character text
% fields, so values must be short: round outward to integers (floor/ceil so
% the rounded range still contains the data -> no clipping) and clamp the
% magnitude so the printed value fits in 8 chars (sign + up to 7 digits).
pmin = floor(min(data, [], 1));
pmax = ceil( max(data, [], 1));
flat = pmax <= pmin;
pmax(flat) = pmin(flat) + 1;
lim = 9999999;
pmin = max(pmin, -lim);
pmax = min(pmax,  lim);

hdr = edfheader("EDF");
hdr.NumDataRecords     = 1;
hdr.DataRecordDuration = seconds(nSamp / fs);
hdr.NumSignals         = nCh;
hdr.SignalLabels       = labs;
hdr.PhysicalDimensions = repelem("uV", nCh);
hdr.PhysicalMin        = pmin;
hdr.PhysicalMax        = pmax;
hdr.DigitalMin         = repmat(-32768, 1, nCh);
hdr.DigitalMax         = repmat( 32767, 1, nCh);

edfw = edfwrite(edf_path, hdr, data, "InputSampleType", "physical"); %#ok<NASGU>

end
