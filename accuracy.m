function [errors,CFL_actual,delt_values] = accuracy(func, n)
% This function calculates the accuracy of an algorithm in relation to an initial state. Note that the CFL number and
% other parameters are set in relation to the total system in CALCULATIONS.m

% func = algorithm of interest
% n = parameter to specify initial conditions

    a = 1;
    tspan = 1;
    L = 1;

    if n == 1
        eqn = @(x) (x >= pi/4 & x <= pi/2);
    else
        k = 1;
        eqn = @(x) (sin(k.*x)).^2;
    end 

    CFL_target = 0.5;
    delt_values = [0.1, 0.05, 0.025, 0.0125];
    errors = zeros(size(delt_values));
    
    % Accuracy Test:
    for i = 1:length(delt_values)
        delt = delt_values(i);
        [U, CFL_actual] = func(a, CFL_target, tspan, delt, L, eqn);
    
        J = size(U, 2);
        x = linspace(0, L, J);
        u_exact = eqn(x);
    
        errors(i) = sqrt(sum((U(end,:) - u_exact).^2) * (L / (J - 1)));
    end



end
