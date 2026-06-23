function clean_labels = decompose_labels(chLabels)
% DECOMPOSE_LABELS  Normalize channel labels (strip leading zeros, etc.).
% Self-contained copy for the realtime_prelim sub-project.

clean_labels = cell(length(chLabels),1);

for ich = 1:length(chLabels)
    if ischar(chLabels)
        label = chLabels;
    else
        label = chLabels{ich};
    end

    if isstring(label)
        label = convertStringsToChars(label);
    end

    label_num_idx = regexp(label,'\d');
    if ~isempty(label_num_idx)
        if ~isscalar(label_num_idx)
            label_num_idx = label_num_idx(1);
        end
        label_non_num = label(1:label_num_idx-1);
        label_num = label(label_num_idx:end);
        if strcmp(label_num(1),'0')
            label_num(1) = [];
        end
        label = [label_non_num,label_num];
    end
    clean_labels{ich} = label;
end

end
