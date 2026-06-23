function [bipolarPairs,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(contacts,contactIndices)
% FIND_BIPOLAR_PAIRS  For each contact, find the next-numbered contact on
% the same electrode (e.g. LA1 -> LA2) to form bipolar pairs.
%
% Outputs:
%   bipolarPairs      : {N x 2} cell of label pairs
%   bipolarIndices    : [N x 2] indices into contactIndices for each pair
%   altBipolarIndices : [numel(contacts) x 1], for each contact the index of
%                       its "+1" partner (NaN if none)

bipolarPairs = {};
bipolarIndices = [];
altBipolarIndices = nan(length(contacts),1);

for i = 1:numel(contacts)
    match = regexp(contacts{i}, '([A-Za-z]+)(\d+)', 'tokens');
    if ~isempty(match)
        letterPart = match{1}{1};
        numberPart = str2double(match{1}{2});
        nextContact = sprintf('%s%d', letterPart, numberPart + 1);
        if any(strcmp(contacts, nextContact))
            bipolarPairs = [bipolarPairs; {contacts{i}, nextContact}]; %#ok<AGROW>
            bipolarIndices = [bipolarIndices; contactIndices(i) contactIndices(find(strcmp(contacts,nextContact),1))]; %#ok<AGROW>
            altBipolarIndices(i) = find(strcmp(contacts,nextContact),1);
        end
    end
end

end
