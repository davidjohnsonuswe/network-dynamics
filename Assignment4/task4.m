%% Task 4 The H1N1 pandemic in Sweden 2009
clear all

N = 10; % Iteration
weeks = 15; % Number of weeks to simulate
n = 934; % Number of nodes
k0 = 10; % Average degree
% Probability from an infected individual to a susceptible one
beta0 = 0.3;
% Probability an infected individual will recover
rho0 = 0.6;
delta_k = 1;
delta_beta = 0.05;
delta_rho = 0.05;

set_k = [k0-delta_k,k0,k0+delta_k];
set_beta = [beta0-delta_beta,beta0,beta0+delta_beta];
set_rho = [rho0-delta_rho,rho0,rho0+delta_rho];

S = 0; % Susceptible 
I = 1; % Infected
R = 2; % Recovered
V = 3; % Vaccinated



% Store optimize values
k_optimize = 0;
beta_optimize = 0;
rho_optimize = 0;

% Vaccine distribution
Vacc = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
% Convert to percentage
Vacc = Vacc/100;
% Number of newly infected individuals each week
I0 = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0]';
count = 0;
% To store RMSE min
RMSE_min = 1000;
while true
    for k_val = set_k
        for beta_val = set_beta
            for rho_val = set_rho
                count = count+1;
                % Generate random graph
                W = random_graph(n,k_val);
                
                % Initialize variables to store results
                newlyInfected = zeros(weeks, N);
                newlyVaccinated = zeros(weeks, N);
                susceptible = zeros(weeks, N);
                infected = zeros(weeks, N);
                recovered = zeros(weeks, N);
                vaccinated = zeros(weeks, N);
                % Run N simulations
                for iter = 1:N
                    X = zeros(n, weeks); % Initialize state matrix
                    infected_nbr = 1;
                    initialInfected = randperm(n, 1); % Select 1 infected nodes at random
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
                            if rand() <= (1 - (1 - beta_val)^length(infected_neighbors))
                                X(susceptibleNodes(i),t+1) = I; % Node becomes infected
                                newlyInfected(t+1,iter) = newlyInfected(t+1,iter)+1;
                            end
                        end
                        % Get indices of infected nodes
                        infectedNodes = find(X(:,t) == I);
                        for i = 1:length(infectedNodes)
                            if rand() <= rho_val
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
                avgNewlyVaccinated = [0;mean(newlyVaccinated,2)];

                RMSE_temp = sqrt((1/15)*sum((avgNewlyInfected-I0).^2));
                if RMSE_temp <= RMSE_min
                    % Stop if the same set of parameters
                    if (k_optimize == k_val)&&(beta_optimize == beta_val)&&(rho_optimize == rho_val)
                        break;
                    end
                    RMSE_min = RMSE_temp;
                    sprintf('Iteration %d, RMSE: %d',count,RMSE_min)
                    k_optimize = k_val;
                    beta_optimize = beta_val;
                    rho_optimize = rho_val;
                    avgNewlyInfected_best = avgNewlyInfected;
                    avgSusceptible_best = avgSusceptible;
                    avgInfected_best = avgInfected;
                    avgRecovered_best = avgRecovered;
                    avgVaccinated_best = avgVaccinated;
                end
            end
        end
    end

% Break if k0 == 3
if k0 == 3
    break
end
k0 = k0-1;
beta0 = beta0-delta_beta;
rho0 = rho0-delta_rho;
set_k = [k0-delta_k,k0,k0+delta_k];
set_beta = [beta0-delta_beta,beta0,beta0+delta_beta];
set_rho = [rho0-delta_rho,rho0,rho0+delta_rho];

end


%% Plot the results
x = 1:weeks;
figure(1)
hold on;
plot(x, avgNewlyInfected_best, 'b-');
plot(x, I0, 'm-');
xlabel('Week');
ylabel('Number of Individuals');
legend('Simulation', 'Real','Location','best');
print('Newly_infected_4.eps','-depsc');

x = 0:weeks;
figure(2)
plot(x, avgSusceptible_best, 'b-');
hold on;
plot(x, avgInfected_best, 'r-');
plot(x, avgRecovered_best, 'm-');
plot(x, avgVaccinated_best, 'k-');
xlabel('Week');
ylabel('Number of Individuals');
legend('Average Susceptible', 'Average Infected', 'Average Recovered', 'Average Vaccinated','Location','best');
%print('Average_4.eps','-depsc');

