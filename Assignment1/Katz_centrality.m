function [sectors, sectors_wrtr] = Katz_centrality(W,name)
    % Remove rows and columns that are all zero
    zeroRows = all(W==0, 2);
    zeroCols = all(W==0, 1);
    W(zeroRows,:) = [];
    W(:,zeroCols) = [];
    name_b = name;
    name_b(zeroRows,:) = [];

    Beta = 0.15;
    N = length(W);
    mu = ones(N,1);

    % mu is 1 for '31 Wholesale & retail trade; repairs'
    mu_wrtr = zeros(N,1);
    idx_wrtr = find(name_b == "31 Wholesale & retail trade; repairs");
    mu_wrtr(idx_wrtr) = 1; 
    
    % Find the eigenvalues and eigenvectors
    % V is the matrix of eigenvectors, and D is the diagonal matrix of eigenvalues.
    [V,D] = eig(W');
    
    % Extracts the diagonal elements of the matrix D
    lambda = diag(D);
    
    [~,i] = sort(lambda,1,'descend');
    t = lambda(i(1));
    % Calculate the Katz centrality
    z = (diag(ones(N,1))-1/lambda(i(1))*(1-Beta)*W')\(Beta*mu);
    % Chosing the top 3 values
    idx = Find_idx(z,3);
    % Return sectors with name
    sectors = [name_b(idx(1)), name_b(idx(2)), name_b(idx(3))];

    % Calculate the Katz centrality with mu_wrtr
    z = (diag(ones(N,1))-1/lambda(i(1))*(1-Beta)*W')\(Beta*mu_wrtr);
    % Chosing the top 3 values
    idx = Find_idx(z,3);
    % Return sectors with name
    sectors_wrtr = [name_b(idx(1)), name_b(idx(2)), name_b(idx(3))];
end