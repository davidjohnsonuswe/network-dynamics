%% Load values
format long
clear all
load('twitter.mat', '-ascii')
load('users.mat', '-ascii')
W = spconvert(twitter);

%% a) Compute iteratively the PageRank and find the five most central nodes.
W(6893,6893) = 0; %Make W square

%Parameters
[n,m]=size(W);
beta=0.15;
mu=ones(n,1);

%Add a self loop to avoid out-degree being zero of first node
out = W*ones(m,1);
for i =1:m 
    if out(i) ==0
        W(i,i) = 1;
    end
end

%Calculate normalized matrix P
w = W*ones(n,1);
D = diag(w);
P = sparse(inv(D)*W);

%Iteratively calculate z with Page-Rank algorithm
iter=500;
z = 0; %initial value

for k=0:iter
    z = z + beta*(1-beta)^k*mu;
    mu = (P'*mu);
end 

%Get 5 most central twitter ids
N=5;
[B, I] = maxk(z,N);

for i=1:N
    idx=I(i);
    twitter_ids(i)=users(idx);
end

%% b) Simulate the discrete-time consensus algorithm with two stubborn nodes, e.g., 
% one with value 0 and one with value 1. You can initiate the rest of the nodes to the 
% neutral value 0.5. Plot how the opinions change over time for a handful of nodes of your choice.

%Let the first two nodes be stubborn
Q=P(3:n, 3:n);
E=P(3:n, 1:2);
%get pagerank
x_=pagerank(n, 500, Q, E, true);

%% c) Investigate how the the choice of nodes with respect to their PageRank change the stationary opinion distribution. Plot the opinion distribution for a couple of choices of the stubborn nodes.

%The third most central and one random node, 9 and 112
stub = [112 3426];
regular = setdiff(1:n, stub);
Q = P(regular, regular);
E = P(regular, stub);
pagerank(n, 500, Q, E, true);

%Two random nodes, 67 and 4552
stub = [67 4552];
regular = setdiff(1:n, stub);
Q = P(regular, regular);
E = P(regular, stub);
pagerank(n, 500, Q, E, true);

%Two random nodes, 167 and 5263
stub = [165 5262];
regular = setdiff(1:n, stub);
Q = P(regular, regular);
E = P(regular, stub);
pagerank(n, 500, Q, E, true);

