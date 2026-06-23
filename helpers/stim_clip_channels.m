function [idx, labs] = stim_clip_channels(chLabels, stimPairLabels, mode)
% STIM_CLIP_CHANNELS  Channels to save for a stim clip.
%
%   [idx, labs] = stim_clip_channels(chLabels, {'LB3','LB4'})            % whole electrode (default)
%   [idx, labs] = stim_clip_channels(chLabels, {'LB3','LB4'}, 'electrode')
%   [idx, labs] = stim_clip_channels(chLabels, {'LB3','LB4'}, 'neighbors')
%
% mode:
%   'electrode' (default) : every contact on the stim electrode (same letter
%                           prefix), in numeric order. This is what the AD
%                           detector needs and keeps clips small.
%   'neighbors'           : just the two stim contacts plus one neighbor on
%                           each side (e.g. {LB2,LB3,LB4,LB5}).
%
% INPUT
%   chLabels       : {nCh x 1} cleaned channel labels (decompose_labels)
%   stimPairLabels : 1x2 cell of the detected bipolar stim pair labels
%
% OUTPUT
%   idx  : indices into chLabels (anatomical order, low->high contact number)
%   labs : the corresponding labels

if nargin < 3 || isempty(mode), mode = 'electrode'; end

% normalize input into a 1x2 cellstr {'LBa','LBb'}
if ischar(stimPairLabels)
    stimPairLabels = strsplit(stimPairLabels, '-');     % 'LB3-LB4' -> {'LB3','LB4'}
elseif isstring(stimPairLabels)
    if isscalar(stimPairLabels)
        stimPairLabels = cellstr(strsplit(stimPairLabels, '-'));
    else
        stimPairLabels = cellstr(stimPairLabels);
    end
end
assert(iscell(stimPairLabels) && numel(stimPairLabels) == 2, ...
    'stimPairLabels must be a 1x2 cell of contact labels, e.g. {''LB3'',''LB4''}.');

% parse letter + number from each stim contact
[letter1, n1] = parse_contact(stimPairLabels{1});
[letter2, n2] = parse_contact(stimPairLabels{2});
assert(strcmp(letter1, letter2), ...
    'Stim pair on different electrodes: %s vs %s', stimPairLabels{1}, stimPairLabels{2});
letter = letter1;

lo = min(n1, n2);
hi = max(n1, n2);

switch lower(mode)
    case 'neighbors'
        % one below, the pair, one above
        cand_nums = [lo-1, lo, hi, hi+1];
        cand_nums = cand_nums(cand_nums >= 1);
    case 'electrode'
        % every contact sharing this electrode's letter prefix
        [allNums, allIdx] = electrode_contacts(chLabels, letter);
        idx  = allIdx;
        labs = chLabels(allIdx);
        % already sorted by contact number below
        [~, order] = sort(allNums);
        idx = idx(order); labs = labs(order);
        return
    otherwise
        error('Unknown mode ''%s'' (use ''electrode'' or ''neighbors'').', mode);
end

idx = [];
labs = {};
for k = 1:numel(cand_nums)
    lbl = sprintf('%s%d', letter, cand_nums(k));
    j = find(strcmp(chLabels, lbl), 1);
    if ~isempty(j)
        idx(end+1,1) = j;       %#ok<AGROW>
        labs{end+1,1} = lbl;    %#ok<AGROW>
    end
end

end


function [nums, idx] = electrode_contacts(chLabels, letter)
% All channels whose letter prefix exactly equals `letter`, with their
% trailing contact numbers.
nums = []; idx = [];
for j = 1:numel(chLabels)
    m = regexp(chLabels{j}, '^([A-Za-z]+)(\d+)$', 'tokens', 'once');
    if ~isempty(m) && strcmp(m{1}, letter)
        idx(end+1,1)  = j;             %#ok<AGROW>
        nums(end+1,1) = str2double(m{2}); %#ok<AGROW>
    end
end
end


function [letter, num] = parse_contact(lbl)
m = regexp(lbl, '([A-Za-z]+)(\d+)', 'tokens', 'once');
assert(~isempty(m), 'Could not parse contact label: %s', lbl);
letter = m{1};
num = str2double(m{2});
end
