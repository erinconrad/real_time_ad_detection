function [pairs, plabels] = ad_candidate_pairs(labels, stim_pair, K)
% AD_CANDIDATE_PAIRS  The K bipolar pairs nearest the stim pair that do NOT
% share a contact with it -- where afterdischarges show up cleanly.
%
%   [pairs, plabels] = ad_candidate_pairs(labels, {'LA3','LA4'}, 2)
%       -> e.g. pairs for LA5-LA6 / ... ; plabels {'LA5-LA6', ...}
%
% This reproduces the heuristic from the annotations: for stim LA5-LA6 it
% returns LA3-LA4 and LA7-LA8 (one each side); for an end pair LA1-LA2 it
% returns LA3-LA4 and LA4-LA5 (both on the available side), simply by taking
% the K nearest non-stim-sharing pairs by contact position.
%
% INPUT
%   labels    : {nCh x 1} electrode contact labels (decompose_labels)
%   stim_pair : 1x2 cell of stim contact labels
%   K         : number of candidate pairs to return (default 2)
%
% OUTPUT
%   pairs   : [m x 2] indices into labels for each bipolar pair (m <= K)
%   plabels : {m x 1} 'LAx-LAy' labels

if nargin < 3 || isempty(K), K = 2; end

nCh = numel(labels);
[bp, bidx] = find_bipolar_pairs(labels(:), 1:nCh);   % all adjacent pairs
nP = size(bp,1);
if nP == 0, pairs = []; plabels = {}; return; end

% stim contact numbers
sn = nan(1,2);
for i = 1:2
    m = regexp(stim_pair{i}, '(\d+)$', 'tokens', 'once');
    if ~isempty(m), sn(i) = str2double(m{1}); end
end
stim_center = mean(sn);

% per-pair: contact numbers, center, whether it shares a stim contact
keep = false(nP,1); center = nan(nP,1);
for r = 1:nP
    n1 = num_of(bp{r,1}); n2 = num_of(bp{r,2});
    center(r) = mean([n1 n2]);
    keep(r) = ~ismember(stim_pair{1}, bp(r,:)) && ~ismember(stim_pair{2}, bp(r,:));
end

idx = find(keep);
[~, ord] = sort(abs(center(idx) - stim_center));   % nearest first
idx = idx(ord);
idx = idx(1:min(K, numel(idx)));

pairs = bidx(idx, :);
plabels = cell(numel(idx),1);
for r = 1:numel(idx)
    plabels{r} = [bp{idx(r),1} '-' bp{idx(r),2}];
end
end

function n = num_of(lbl)
m = regexp(lbl, '(\d+)$', 'tokens', 'once');
n = str2double(m{1});
end
