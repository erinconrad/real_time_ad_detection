function exclude = find_exclude_chs(chLabels)
% FIND_EXCLUDE_CHS  Logical mask of scalp/reference/EKG channels to ignore.
excluded = {'C3','C4','CZ','EKG1','EKG2','FZ','LOC','ROC'};
exclude = ismember(chLabels,excluded);
end
