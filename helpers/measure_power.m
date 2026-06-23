function p = measure_power(data,freq_band,fs)
% MEASURE_POWER  Power of a data chunk (samples x channels).
%   measure_power(data)            -> sum of squares (line-length-like energy)
%   measure_power(data,band,fs)    -> bandpower in the given [lo hi] band
%
% NaNs are replaced by the per-channel mean before computing.

meanValues = mean(data,1,'omitnan');
nanIndices = isnan(data);
for col = 1:size(data,2)
    data(nanIndices(:,col),col) = meanValues(col);
end

if nargin < 2 || isempty(freq_band)
    p = sum(data.^2,1);
else
    data(isnan(data)) = 0;
    p = bandpower(data,fs,freq_band);
end

end
