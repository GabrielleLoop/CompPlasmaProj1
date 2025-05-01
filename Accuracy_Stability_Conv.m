% === Parameters ===
clear variables
close all
clc




[errors1,~,delt_values1] = accuracy(@FTBS, 1);
[errors2,~,delt_values2] = accuracy(@FTCS, 1);
[errors3,~,delt_values3] = accuracy(@LF, 1);
[errors4,~,delt_values4] = accuracy(@LW, 1);
[errors5,CFL_actual1,delt_values5] = accuracy(@HighOrder, 1);

% === Plot Accuracy ===
figure(1);
loglog(delt_values1, errors1, '-o', 'LineWidth', 1);
hold on
loglog(delt_values2, errors2, '-o', 'LineWidth', 1);
loglog(delt_values3, errors3, '-o', 'LineWidth', 1);
loglog(delt_values4, errors4, '-o', 'LineWidth', 1);
loglog(delt_values5, errors5, '-o', 'LineWidth', 1);
xlabel('\Deltat'); ylabel('L2 Error');
title('Accuracy Test with Square Function, CFL = 0.4875');
legend('FTBS','FTCS','LF','LW','High Order')
grid on;
hold off

[errors1,~,delt_values1] = accuracy(@FTBS, 2);
[errors2,~,delt_values2] = accuracy(@FTCS, 2);
[errors3,~,delt_values3] = accuracy(@LF, 2);
[errors4,~,delt_values4] = accuracy(@LW, 2);
[errors5,CFL_actual2,delt_values5] = accuracy(@HighOrder, 2);

% === Plot Accuracy ===
figure(2);
loglog(delt_values1, errors1, '-o', 'LineWidth', 1);
hold on
loglog(delt_values2, errors2, '-o', 'LineWidth', 1);
loglog(delt_values3, errors3, '-o', 'LineWidth', 1);
loglog(delt_values4, errors4, '-o', 'LineWidth', 1);
loglog(delt_values5, errors5, '-o', 'LineWidth', 1);
xlabel('\Deltat'); ylabel('L2 Error');
title('Accuracy Test with Sine Squared Function, CFL = 0.4875');
legend('FTBS','FTCS','LF','LW','High Order')
grid on;
hold off




%%
clc
clear vars;
close all;

n=1;
func = FTBS;

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





% === Plot Accuracy ===
figure(1);
loglog(delt_values, errors, '-o', 'LineWidth', 2);
xlabel('\Deltat'); ylabel('L2 Error');
title('FTBS Accuracy Test with External Function');
grid on;



%%













% === Estimate Order ===
p = polyfit(log(delt_values), log(errors), 1);
fprintf('Estimated order of accuracy: %.2f\n', p(1));

% === Stability Test ===
CFL_values = [0.5, 1.0, 1.2];
dx = 0.01;
delt = CFL_values(1) * dx / a;

figure(2); 
hold on
for k = 1:length(CFL_values)
    CFL_k = CFL_values(k);
    delt_k = CFL_k * dx / a;

    [U, ~] = FTBS(a, CFL_k, tspan, delt_k, L, initial_fn);
    x = linspace(0, L, size(U, 2));
    plot(x, U(end,:), 'DisplayName', ['CFL = ' num2str(CFL_k)], 'LineWidth', 1.5);
end
xlabel('x'); ylabel('u');
title('FTBS Stability Test with External Function');
legend show; grid on;