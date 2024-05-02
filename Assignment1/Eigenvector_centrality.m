function sectors = Eigenvector_centrality(W,name)
    % Remove rows and columns that are all zero
    zeroRows = all(W==0, 2);
    zeroCols = all(W==0, 1);
    W(zeroRows,:) = [];
    W(:,zeroCols) = [];
    name_b = name;
    name_b(zeroRows,:) = [];
    
    % Find the eigenvalues and eigenvectors
    % V is the matrix of eigenvectors, and D is the diagonal matrix of eigenvalues.
    [V,D] = eig(W');
    
    % Extracts the diagonal elements of the matrix D
    lambda = diag(D);
    
    [~,i] = sort(lambda,1,'descend');
    % The eigenvector centrality
    z = abs(V(:,i(1)));
    % Chosing the top 3 values in the eigenvector
    idx = Find_idx(z,3);
    % Return sectors with name
    sectors = [name_b(idx(1)), name_b(idx(2)), name_b(idx(3))];
end