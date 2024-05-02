%% 2a) Simulate the learning dynamics
clear all

% Define the adjacency matrix
W = [0 1 0 0 0 0 0 0 0 0;
       1 0 1 0 0 0 0 0 0 0;
       0 1 0 1 0 0 0 0 0 0;
       0 0 1 0 1 0 0 0 0 0;
       0 0 0 1 0 1 0 0 0 0;
       0 0 0 0 1 0 1 0 0 0;
       0 0 0 0 0 1 0 1 0 0;
       0 0 0 0 0 0 1 0 1 0;
       0 0 0 0 0 0 0 1 0 1;
       0 0 0 0 0 0 0 0 1 0];

nbr_nodes = 10;
% 1 is red, 2 is green
node_colors = [1 2];
% Initialize the potential
potential = 1;
% Potential array to store data
potential_arr = zeros(1,1);
% Initialize every node to red
X = ones(10,1);
% Node coordinates
coordinate = [0:10; zeros(1, 11)]';
t = 0;

% Plot
figure(1)
gplot(W,coordinate,'-k');
hold on
colors = ['r', 'g'];
for i = 1:nbr_nodes
    scatter(coordinate(i,1), coordinate(i,2), 100, colors(X(i)), 'filled', 'markeredgecolor', 'k');
end
% Remove the x-axis and y-axis tick marks
set(gca, 'xtick', [], 'ytick', [], 'Color', 'none')
% Get the current position of the figure
figure_position = get(gcf, 'Position');
% Reduce the height to one-third of the default value
figure_position(4) = figure_position(4) / 3;
% Modify the size of the figure
set(gcf, 'Position', figure_position);
%print('Initialize_all_nodes_to_red.eps','-depsc');

% Loop until the potential is zero
while potential ~= 0
    t = t+1;
    eta = t/100;
    potential = 0;
    % Choose a random node
    node = randi(nbr_nodes);

    % Calculate the probability distribution
    prob_distribution = calc_probdistribution(W,node,node_colors,eta,X);
    color_value = randsample(node_colors,1,true,prob_distribution);
    X(node) = color_value;
   
    
    % Calculate the potential then add to the potential array
    for i=1:nbr_nodes
        for j = find(W(i,:) ~= 0)
            potential = potential + (1/2)*cost_function(X(i),X(j));
        end
    end
    potential_arr(t) = potential;
end

% Plot
figure(2)
gplot(W,coordinate,'-k');
hold on
colors = {'r', 'g'};
for i = 1:nbr_nodes
    scatter(coordinate(i,1), coordinate(i,2), 100, colors{X(i)}, 'filled', 'markeredgecolor', 'k');
end
% Remove the x-axis and y-axis tick marks
set(gca, 'xtick', [], 'ytick', [], 'Color', 'none')
% Get the current position of the figure
figure_position = get(gcf, 'Position');
% Reduce the height to one-third of the default value
figure_position(4) = figure_position(4) / 3;
% Modify the size of the figure
set(gcf, 'Position', figure_position);
%print('Assign_colors_to_nodes.eps','-depsc');

figure(3)
plot(potential_arr)
xlabel('Time')
ylabel('Potential function')
%print('Potential_function.eps','-depsc');

%% Helper functions
function c = cost_function(s,X_t)
    c = double(s == X_t);
end

function prob_distribution = calc_probdistribution(W,node,node_colors,eta,X)
    % Probability distribution of nodes
    prob_distribution = zeros(length(node_colors),1);
    for color = node_colors
        sum_c = 0;
        % Calculate the sum of c(color,X_j(t))
        for j = find(W(node,:)~=0)
            sum_c = sum_c + cost_function(color,X(j));
        end
        prob_distribution(color) = exp(-eta*sum_c);
    end
    % Calculate the probability distribution
    prob_distribution = prob_distribution/sum(prob_distribution);
end