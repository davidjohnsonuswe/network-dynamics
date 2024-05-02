%% Load data
load -ascii flow.mat
load -ascii traffic.mat
load -ascii capacities.mat
load -ascii traveltime.mat

% The node-link incidence matrix B
B = traffic;
% Find the number of nodes and links
[num_nodes, num_links] = size(B);
% Maximum flow capacity
C = capacities;
% Minimum traveling times
l = traveltime;

%% a) Find the shortest path between node 1 and 17

% Initialize source and target matrices with zeros
source = zeros(num_links, 1);
target = zeros(num_links, 1);

% Loop through each link in B
for i = 1:num_links
    % Find the column indices of the non-zero elements
    col_indices = find(B(:, i) ~= 0);
    
    % Set the source and target nodes
    source(i) = col_indices(1);
    target(i) = col_indices(2);
end

G = digraph(source,target,l);
[P,d] = shortestpath(G,1,17);

%% b) Find the maximum flow between node 1 and 17.
Gb = digraph(source,target,C);
mf = maxflow(Gb,1,17);

%% c) compute the external inflow or outflow at each node
in_out_flow = B*flow;

%% d) Find the social optimum f∗ with respect to the delays on the different links de(fe)

% Exogenous flow
nu = zeros(num_nodes, 1);
nu(1) = in_out_flow(1);
nu(num_nodes) = -in_out_flow(1);

% Compute the social optimum
cvx_begin
    variable f(num_links);
    minimize sum(l.*C.*inv_pos(1-f.*inv_pos(C))-l.*C);
    subject to
        B*f == nu;
        0 <= f <= C;
cvx_end
% Print the flow vector
f

%% e) Find the Wardrop equilibrium f(0)

% Compute the Wardrop equilibrium
cvx_begin
    variable f0(num_links);
    minimize sum(-l.*C.*log(1-f0.*inv_pos(C)));
    subject to
        B*f0 == nu;
        0 <= f0 <= C;
cvx_end
% Print the flow vector
f0

%% f) Compute the new Wardrop equilibrium f(ω)

% fstar is the flow at the system optimum
fstar = f;
% Calculate w
w = fstar.*((l.*C)./((C-fstar).^2));

% Compute the Wardrop equilibrium
cvx_begin
    variable fw(num_links);
    minimize sum(-l.*C.*log(1-fw.*inv_pos(C)) + fw.*w);
    subject to
        B*fw == nu;
        0 <= fw <= C;
cvx_end
% Print the flow vector
fw
f_diff = fw - f;

%% g) Compute the newWardrop equilibrium with the constructed tolls f(ω∗)

% Compute the system optimum with the new cost function
cvx_begin
    variable f_n(num_links);
    minimize sum(l.*C.*inv_pos(1-f_n.*inv_pos(C))-l.*C-f_n.*l);
    subject to
        B*f_n == nu;
        0 <= f_n <= C;
cvx_end
% Print the flow vector
f_n

% Construc w_star
f_star_g = f_n;
w_star = 1./quad_over_lin(C-f_star_g,f_star_g.*l.*C,num_links) - l;

% Compute the Wardrop equilibrium
cvx_begin
    variable fw_star(num_links);
    minimize sum(-l.*C.*log(1-fw_star.*inv_pos(C)) + fw_star.*w_star);
    subject to
        B*fw_star == nu;
        0 <= fw_star <= C;
cvx_end
% Print the flow vector
fw_star
f_diff_g = fw_star - f_n;
