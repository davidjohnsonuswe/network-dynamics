%% Task 2 Simulate a pandemic without vaccination
clear all

% Probability from an infected individual to a susceptible one
beta = 0.3;
% Probability an infected individual will recover
rho = 0.7;
N = 100; % Iteration
weeks = 15; % Number of weeks to simulate
k = 6; % Average degree
n = 500; % Number of nodes

% The weight matrix
W = random_graph(n,k);

S = 0; % Susceptible 
I = 1; % Infected
R = 2; % Recovered

% Initialize variables to store results
newlyInfected = zeros(weeks, N);
susceptible = zeros(weeks, N);
infected = zeros(weeks, N);
recovered = zeros(weeks, N);

% Run N simulations
for iter = 1:N
    X = zeros(n, weeks); % Initialize state matrix
    infected_nbr = 10;
    initialInfected = randperm(n, 10); % Select 10 infected nodes at random
    X(initialInfected,1) = I; % Set the initial infected nodes
    newlyInfected(1,iter) = infected_nbr; % Initial newly infected
    susceptible(1,iter) = n-infected_nbr;
    infected(1,iter) = infected_nbr;
    
    % Simulate the epidemic for each week
    for t = 1:weeks-1
        % Simulate disease propagation
        susceptibleNodes = find(X(:,t) == S); % Get indices of susceptible nodes
        for i = 1:length(susceptibleNodes)
            neighbors = find(W(susceptibleNodes(i), :)); % Get indices of neighbors
            infected_neighbors = neighbors((X(neighbors,t)==I));
            if rand() <= (1 - (1 - beta)^length(infected_neighbors))
                X(susceptibleNodes(i),t+1) = I; % Node becomes infected
                newlyInfected(t+1,iter) = newlyInfected(t+1,iter)+1;
            end
        end
        % Get indices of infected nodes
        infectedNodes = find(X(:,t) == I);
        for i = 1:length(infectedNodes)
            if rand() <= rho
                X(infectedNodes(i),t+1) = R; % Node i recovers
            else
                X(infectedNodes(i),t+1) = I;
            end
        end
        % Get indices of recovered nodes
        recoveredNodes = find(X(:,t) == R);
        for i = 1:length(recoveredNodes)
            X(recoveredNodes(i),t+1) = R; % Node i recovers
        end

        % Count the number of susceptible, infected, and recovered individuals
        susceptible(t+1,iter) = length(find(X(:,t+1) == S));
        infected(t+1,iter) = length(find(X(:,t+1) == I));
        recovered(t+1,iter) = length(find(X(:,t+1) == R));
    end
    
end

% Calculate the average
avgNewlyInfected = [mean(newlyInfected,2)];
avgSusceptible = [n;mean(susceptible,2)];
avgInfected = [0;mean(infected,2)];
avgRecovered = [0;mean(recovered,2)];

% Plot the results
x = 1:weeks;
figure(1)
plot(x, avgNewlyInfected, 'b-');
xlabel('Week')
ylabel('Number of newly infected individuals')
print('Newly_infected_2.eps','-depsc');

x = 0:weeks;
figure(2)
plot(x, avgSusceptible, 'b-');
hold on;
plot(x, avgInfected, 'r-');
plot(x, avgRecovered, 'm-');
xlabel('Week');
ylabel('Number of Individuals');
legend('Average Susceptible', 'Average Infected', 'Average Recovered','Location','best');
%print('Average_2.eps','-depsc');

