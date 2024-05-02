function W = random_graph(n,k)

% Generate random graph according to the preferential attachment model

% Number of nodes
k0 = k+1;
% Adjacency matrix
W = ones(k0);
W = triu(W, 1) + triu(W, 1)';

% Preferential attachment
for t = k0+1:n
    % The number of links we wish to add from the new node
    if mod(t,2)==0
        c = floor(k/2);
    else
        c = ceil(k/2);
    end
    w = sum(W,2); % Degree vector of graph
    P = w./sum(w); % Probability vector for adding links
    for j=1:c
        % randsample() will draw one random sample from the population
        % 1:(k0+1). The probability that it draws a certain individual is
        % specified in P. Note that randsample() use replacement, so we have
        % to remove the drawn individual from the population.
        neighbour = randsample(1:k0,1,true,full(P));
        P(neighbour) = 0;
        W(k0+1,neighbour) = 1;
        W(neighbour,k0+1) = 1;
    end
    % Increase number of nodes
    k0 = k0 + 1;
end
end