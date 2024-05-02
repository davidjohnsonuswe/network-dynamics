%% 2b) Assign wifi-channels to routers
clear all

load -ascii data/wifi.mat
load -ascii data/coord.mat

% Adjacency matrix
W_b = wifi;
% Node coordinates
coordinate_b = coord;
nbr_nodes_b = 100;
% 1-red, 2-green, 3-blue, 4-yellow, 5-magenta, 6-cyan, 7-white, 8-black
node_colors_b = [1 2 3 4 5 6 7 8];
% Initialize the potential
potential_b = 1;
% Potential array to store data
potential_arr_b = zeros(1,1);
% Initialize every node to red
X_b = ones(100,1);
t_b = 0;

% Loop until the potential is zero
for k = 1:800
    t_b = t_b+1;
    eta = t_b/300;
    potential_b = 0;
    % Choose a random node
    node_b = randi(nbr_nodes_b);

    % Calculate the probability distribution
    prob_distribution_b = calc_probdistribution_b(W_b,node_b,node_colors_b,eta,X_b);
    color_value_b = randsample(node_colors_b,1,true,prob_distribution_b);
    X_b(node_b) = color_value_b;
   
    
    % Calculate the potential then add to the potential array
    for i=1:nbr_nodes_b
        for j = find(W_b(i,:) ~= 0)
            potential_b = potential_b + (1/2)*cost_function_b(X_b(i),X_b(j));
        end
    end
    potential_arr_b(t_b) = potential_b;
end

% Plot
figure(1)
gplot(W_b,coordinate_b,'-k');
hold on
colors_b = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'];

for i = 1:nbr_nodes_b
    colorIndex = X_b(i, 1);
    scatter(coordinate_b(i,1), coordinate_b(i,2), 100, colors_b(colorIndex), 'filled', 'markeredgecolor', 'k');
end
% Remove the x-axis and y-axis tick marks
set(gca, 'xtick', [], 'ytick', [], 'Color', 'none')
print('Assign_colors_to_nodes_b.eps','-depsc');

figure(2)
plot(potential_arr_b)
xlabel('Time')
ylabel('Potential function')
print('Potential_function_b.eps','-depsc');

%% Helper functions
function c = cost_function_b(s,X_t)
    if s == X_t
        c = 2;
    elseif abs(X_t-s) == 1
        c = 1;
    else
        c = 0;
    end
end

function prob_distribution_b = calc_probdistribution_b(W,node,node_colors,eta,X)
    % Probability distribution of nodes
    prob_distribution = zeros(length(node_colors),1);
    for color = node_colors
        sum_c = 0;
        % Calculate the sum of c(color,X_j(t))
        for j = find(W(node,:)~=0)
            sum_c = sum_c + cost_function_b(color,X(j));
        end
        prob_distribution(color) = exp(-eta*sum_c);
    end
    % Calculate the probability distribution
    prob_distribution_b = prob_distribution/sum(prob_distribution);
end