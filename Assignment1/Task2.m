%% Load data
load('twitter.mat', '-ascii')
load('users.mat', '-ascii')
W = spconvert(twitter);

%% 2a) Compute iteratively the PageRank and find the five most central nodes

% Make W a square matrix
W(6893,6893) = 0;

n = length(W);
beta = 0.15;
mu = ones(n,1);

% Check the out-degrees of all nodes
out_degree = W*ones(n,1);
for i = 1:n 
    if out_degree(i) == 0
        W(i,i) = 1;
    end
end


% Calculate the out-degrees
w = W*ones(n,1);
% Diagonal matrix
D = diag(w);
% Calculate the normalized matrix P
P = sparse(inv(D)*W);

z = 0;

% Calculate the page rank centrality vector z.
for k=0:100
    z = z + beta*(1-beta)^k*mu;
    mu = (P'*mu);
end

% tol = 1e-6;  % convergence threshold
% diff = inf;  % initialize difference
% k = 0;       % initialize iteration counter
% 
% while diff > tol
%     z_new = beta * P' * z + (1 - beta) * mu;  % update PageRank vector
%     diff = norm(z_new - z, 1);                % compute difference
%     z = z_new;                                % update PageRank vector
%     k = k + 1;                                % increment iteration counter
% end

% Chosing the top 5 values in z
idx = Find_idx(z,5);
% Find the index of the min value. Use in 2c
min_idx = find(z == min(z), 1, 'first');
% Return sectors with name
sectors = [users(idx(1)), users(idx(2)), users(idx(3)), users(idx(4)), users(idx(5))];

disp("The five most central nodes")
disp(idx)

%% 2b) Simulate the discrete-time consensus algorithm with two stubborn nodes

% Select the first 2 nodes as the stubborn nodes
Q = P(3:n,3:n);
E = P(3:n,1:2);

% Node 1 has value 1 and node 2 has value 0
u = [1;0];
% Initiate the rest with 0.5
x_u = 1/2*ones(length(Q),1);
x = [x_u;u];

% Choose some nodes to calculate values.
node5 = zeros(500,1);
node900 = zeros(500,1);
node5000 = zeros(500,1);
for i = 1:500
    x_u = Q*x_u + E*u;
    x = [x_u;u];
    node5(i) = x_u(5);
    node900(i) = x_u(900);
    node5000(i) = x_u(5000);
end

figure(1)
plot(1:500,node5)
hold on
plot(1:500,node900)
hold on
plot(1:500,node5000)
xlabel('Number of steps');
ylabel('Opinion value')
legend('node5', 'node900','node5000','Location','best');
%print('Opinions_over_time.eps','-depsc');


figure(2)
histogram(x)
xlabel('Opinion value');
ylabel('Number of nodes')
%print('Histogram_opinions.eps','-depsc');

%% 2c) Investigate how the choice of nodes with respect to their PageRank
% change the stationary opinion distribution

% Choose the third most central and one random node 100
stb_nodes = [idx(3) 100];
regular_nodes = setdiff(1:n, stb_nodes);
Q = P(regular_nodes, regular_nodes);
E = P(regular_nodes, stb_nodes);
x = Distributed_averaging(Q,E);

figure(3)
histogram(x)
xlabel('Opinion value');
ylabel('Number of nodes')
%print('2c1.eps','-depsc');

% Choose the least central and one random node 1000
stb_nodes = [min_idx 1000];
regular_nodes = setdiff(1:n, stb_nodes);
Q = P(regular_nodes, regular_nodes);
E = P(regular_nodes, stb_nodes);
x = Distributed_averaging(Q,E);

figure(4)
histogram(x)
xlabel('Opinion value');
ylabel('Number of nodes')
%print('2c2.eps','-depsc');

% Choose the 2 random nodes 2000 and 3000
stb_nodes = [2000 3000];
regular_nodes = setdiff(1:n, stb_nodes);
Q = P(regular_nodes, regular_nodes);
E = P(regular_nodes, stb_nodes);
x = Distributed_averaging(Q,E);

figure(5)
histogram(x)
xlabel('Opinion value');
ylabel('Number of nodes')
%print('2c3.eps','-depsc');
