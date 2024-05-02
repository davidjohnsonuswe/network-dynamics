%% Task 3 Simulate a pandemic with vaccination
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
V = 3; % Vaccinated

% Initialize variables to store results
newlyInfected = zeros(weeks, N);
newlyVaccinated = zeros(weeks, N);
susceptible = zeros(weeks, N);
infected = zeros(weeks, N);
recovered = zeros(weeks, N);
vaccinated = zeros(weeks, N);

% Vaccine distribution
Vacc = [0, 5, 15, 25, 35, 45, 55, 60, 60, 60, 60, 60, 60, 60, 60];
% Convert to percentage
Vacc = Vacc/100;

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
        % Get indices of not vaccinated and vaccinated nodes
        not_vaccinatedNodes = find(~(X(:,t) == V));
        vaccinatedNodes = find((X(:,t) == V));
        % Find number of individuals to vaccinate
        nrb_vacc = uint8((Vacc(t+1)-Vacc(t))*n);
        % Select individuals to vaccinate at random from not_vaccinated
        vacc_idx = randperm(length(not_vaccinatedNodes), nrb_vacc);

        for i = 1:nrb_vacc
            X(vacc_idx(i),t+1) = V; % Node is vaccinated
            newlyVaccinated(t+1,iter) = newlyVaccinated(t+1,iter) + 1;
        end
        % Update vaccinated list for next week
        for i = 1:length(vaccinatedNodes)
            X(vaccinatedNodes(i),t+1) = V;
        end

        susceptibleNodes = find(X(:,t) == S); % Get indices of susceptible nodes
        for i = 1:length(susceptibleNodes)
            neighbors = find(W(susceptibleNodes(i),:)); % Get indices of neighbors
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

        % Count the number of susceptible, infected, vaccinated and recovered individuals
        susceptible(t+1,iter) = length(find(X(:,t+1) == S));
        infected(t+1,iter) = length(find(X(:,t+1) == I));
        recovered(t+1,iter) = length(find(X(:,t+1) == R));
        vaccinated(t+1,iter) = length(find(X(:,t+1) == V));
    end
    
end

% Calculate the average
avgNewlyInfected = [mean(newlyInfected,2)];
avgSusceptible = [n;mean(susceptible,2)];
avgInfected = [0;mean(infected,2)];
avgRecovered = [0;mean(recovered,2)];
avgVaccinated = [0;mean(vaccinated,2)];
avgNewlyVaccinated = [mean(newlyVaccinated,2)];

% Plot the results
x = 1:weeks;
figure(1)
hold on;
plot(x, avgNewlyInfected, 'b-');
plot(x, avgNewlyVaccinated, 'm-');
xlabel('Week');
ylabel('Number of Individuals');
legend('Newly infected', 'Newly vaccinated','Location','best');
print('Newly_infected_vaccinated.eps','-depsc');

x = 0:weeks;
figure(2)
plot(x, avgSusceptible, 'b-');
hold on;
plot(x, avgInfected, 'r-');
plot(x, avgRecovered, 'm-');
plot(x, avgVaccinated, 'k-');
xlabel('Week');
ylabel('Number of Individuals');
legend('Average Susceptible', 'Average Infected', 'Average Recovered', 'Average Vaccinated','Location','best');
%print('Average_3.eps','-depsc');

