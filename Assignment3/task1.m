%% Calculate a transition probability matrix
clear all

% Transition rate matrix
Lambda = [0 2/5 1/5 0 0;
          0 0 3/4 1/4 0;
          1/2 0 0 1/2 0;
          0 0 1/3 0 2/3;
          0 1/3 0 1/3 0];

[M,N] = size(Lambda);

% The rate of the distribution
w = Lambda*ones(M,1);

% Calculate the normalized weight matrix
D = diag(w);
P = D\Lambda;

w_star = max(w);

% Calculate P_bar_ij for i != j
P_bar = zeros(M,N);
for i=1:M
    for j=1:N
       if i~=j
           P_bar(i,j) = Lambda(i,j)/w_star;
       end
    end
end

% Calculate P_bar_ii for i = j
for i = 1:M
    P_bar(i, i) = 1 - sum(P_bar(i, setdiff(1:M, i)));
end

%% a) Average time to return to node a

nodes = [1 2 3 4 5]; % nodes o a b c d
start_node = 2; % node a
iter = 1e5;
% To add the time that the particle returns to node a
time_return = zeros(iter,1);
next_pos = start_node;
t = 0;

for i=1:iter
    
    % Uniform distribution
    u = rand;
    r = w_star;
    t_next = -log(u)/r;
    t = t+t_next;

    % Select the next position based on P_bar
    next_pos = randsample(nodes,1,true,P_bar(next_pos,:));
    
    % Checking if the particle returns to node a
    if next_pos==start_node
        time_return(i) = t;
        t = 0;
    end
end
% Calculate the average time
average_time = mean(time_return(time_return ~= 0));

%% b) Theoretical return-time

[V,D] = eig(P_bar');
% Find the eigenvector corresponding to the eigenvalue 1
pi_bar = V(:, find(abs(diag(D) - 1) < 1e-10));
% Normalize the invariant probability vector
pi_bar = pi_bar/sum(pi_bar);

% Calculate the return time
return_time = 1/(w(2)*pi_bar(2));

%% c) Average time to move from node o to node d

nodes = [1 2 3 4 5]; % nodes o a b c d
start_node = 1; % node o
end_node = 5; % node d
iter = 1e5;
% To add the time that the particle returns to node a
move_time = zeros(iter,1);
next_pos = start_node;
t = 0;

for i=1:iter
    
    % Uniform distribution
    u = rand;
    r = w_star;
    t_next = -log(u)/r;
    t = t+t_next;

    % Select the next position based on P_bar
    next_pos = randsample(nodes, 1, true, P_bar(next_pos,:));
    
    % Checking if the particle reaches node d
    if next_pos==end_node
        move_time(i) = t;
        next_pos = start_node;
        t = 0;
    end
end
% Calculate the average time
average_move_time = mean(move_time(move_time ~= 0));


%% d) Theoretical hitting time from o to d

node_list = [1,2,3,4];
% Remove d(5) from P_bar
P_bar_new = P_bar(node_list,node_list);

A = (eye(length(node_list))-P_bar_new);
B = ones(length(node_list),1);
X = A\B;
% Show the hitting time from node o to node d.
X(1)