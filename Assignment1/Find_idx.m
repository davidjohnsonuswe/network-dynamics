function idx_sector = Find_idx(w,k)
    % sort w in descending order and get sorted values
    sortedValues = sort(w, 'descend');
    idx_sector = zeros(1,k);
    for i = 1:k
        % get the original index of the top k values
        idx_sector(i) = find(w == sortedValues(i), 1, 'first');
    end
end