function x = Distributed_averaging(Q,E)
    % Values of stubborn nodes
    u = [1;0];
    % Initiate the rest with 0.5
    x_u = 1/2*ones(length(Q),1);
    x = [x_u;u];

    for i = 1:500
        x_u = Q*x_u + E*u;
        x = [x_u;u];
    end
end