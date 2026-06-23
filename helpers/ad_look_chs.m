function chs_to_look = ad_look_chs(chLabels,altBipolarIndices,stimChs,file_name)
% AD_LOOK_CHS  Logical mask of channels eligible for AD detection for a
% given stim pair: contacts on the same electrode as the stim channels,
% excluding the stim contacts themselves and (except for HUP260) excluding
% contacts numbered >= 11 (typically outside brain / noisy).

if nargin < 4, file_name = ''; end

all_chs = (1:length(chLabels))';
chs_to_look = false(length(chLabels),1);

% Channels that ARE the stim pair (either member, via bipolar partner)
is_stim_ch = (ismember(all_chs,stimChs) | ismember(altBipolarIndices,stimChs));

% Electrode letter of the stim channel
stimChLabel = chLabels{stimChs(1)};
match = regexp(stimChLabel, '([A-Za-z]+)(\d+)', 'tokens');
letterPart = match{1}{1};

checkLabel = @(x) (length(x) >= length(letterPart)) && ...
    strcmp(letterPart, x(1:length(letterPart)));
same_elec = cellfun(checkLabel, chLabels);

% Trailing contact number for each channel
numericParts = cellfun(@(s) regexp(s, '\d+$', 'match'), chLabels, 'UniformOutput', false);
numericArray = nan(size(chLabels));
notEmptyIdx = ~cellfun(@isempty, numericParts);
numericArray(notEmptyIdx) = cellfun(@(x) str2double(x{1}), numericParts(notEmptyIdx));
less_than_eleven = numericArray < 11;

if strcmp(file_name,'HUP260_phaseII')
    chs_to_look(same_elec & ~is_stim_ch) = true;
else
    chs_to_look(same_elec & ~is_stim_ch & less_than_eleven) = true;
end

% Drop excluded (scalp/ekg) channels
chs_to_look(find_exclude_chs(chLabels)) = false;

end
