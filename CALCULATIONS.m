%% SINE-SQUARED WAVE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section references all the algorithm functions and compiles them to
% overlay their profiles for different CFL numbers and time spans. 

% u(x) = sin^2(kx), k=1
clc
clear variables
close all


% Defining equation:
k = 1;
eqn = @(x) (sin(k.*x)).^2;

% Defining solving parameters:
a = 1;
delt = pi/650;
L = pi;

CFLA = 0.50;
CFLB = 0.75;
CFLC = 0.98;
tspana = pi; % (period*pi)
Na = tspana/delt;
tspanb = 10*pi;
Nb = tspanb/delt;


[U1a, ~] = FTBS(a, CFLA, tspana, delt, L ,eqn);
[U2a, ~] = FTCS(a, CFLA, tspana, delt, L ,eqn);
[U3a, ~] = LF(a, CFLA, tspana, delt, L ,eqn);
[U4a, ~] = LW(a, CFLA, tspana, delt, L ,eqn);
[U5a, ~] = HighOrder(a, CFLA, tspana, delt, L ,eqn);

[U1b, ~] = FTBS(a, CFLA, tspanb, delt, L ,eqn);
[U2b, ~] = FTCS(a, CFLA, tspanb, delt, L ,eqn);
[U3b, ~] = LF(a, CFLA, tspanb, delt, L ,eqn);
[U4b, CFLA] = LW(a, CFLA, tspanb, delt, L ,eqn);
[U5b, ~] = HighOrder(a, CFLA, tspanb, delt, L ,eqn);
plotansCFLA = plotUsin2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
sgtitle('CFL = 0.50')
hold off

[U1a, ~] = FTBS(a, CFLB, tspana, delt, L ,eqn);
[U2a, ~] = FTCS(a, CFLB, tspana, delt, L ,eqn);
[U3a, ~] = LF(a, CFLB, tspana, delt, L ,eqn);
[U4a, ~] = LW(a, CFLB, tspana, delt, L ,eqn);
[U5a, ~] = HighOrder(a, CFLB, tspana, delt, L ,eqn);

[U1b, ~] = FTBS(a, CFLB, tspanb, delt, L ,eqn);
[U2b, ~] = FTCS(a, CFLB, tspanb, delt, L ,eqn);
[U3b, ~] = LF(a, CFLB, tspanb, delt, L ,eqn);
[U4b, CFLB] = LW(a, CFLB, tspanb, delt, L ,eqn);
[U5b, ~] = HighOrder(a, CFLB, tspanb, delt, L ,eqn);
plotansCFLB = plotUsin2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
sgtitle('CFL = 0.75')
hold off

[U1a, ~] = FTBS(a, CFLC, tspana, delt, L ,eqn);
[U2a, ~] = FTCS(a, CFLC, tspana, delt, L ,eqn);
[U3a, ~] = LF(a, CFLC, tspana, delt, L ,eqn);
[U4a, ~] = LW(a, CFLC, tspana, delt, L ,eqn);
[U5a, ~] = HighOrder(a, CFLC, tspana, delt, L ,eqn);

[U1b, ~] = FTBS(a, CFLC, tspanb, delt, L ,eqn);
[U2b, ~] = FTCS(a, CFLC, tspanb, delt, L ,eqn);
[U3b, ~] = LF(a, CFLC, tspanb, delt, L ,eqn);
[U4b, CFLC] = LW(a, CFLC, tspanb, delt, L ,eqn);
[U5b, ~] = HighOrder(a, CFLC, tspanb, delt, L ,eqn);
plotansCFLC = plotUsin2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
sgtitle('CFL = 0.98')
hold off


%% SQUARE WAVE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section does the same procedure as the one above, but the plotting
% of these profiles are different functions so that the "real" solution
% reflects the correct equation.

% u(x) = if pi/4
clc
clear variables
close all


% Defining equation:
eqn = @(x) (x >= pi/4 & x <= pi/2);

% Defining solving parameters:
a = 1;
delt = pi/650;
L = pi;

CFLA = 0.50;
CFLB = 0.75;
CFLC = 0.98;
tspana = pi; % (period*pi)
Na = tspana/delt;
tspanb = 10*pi;
Nb = tspanb/delt;


[U1a, ~] = FTBS(a, CFLA, tspana, delt, L ,eqn);
[U2a, ~] = FTCS(a, CFLA, tspana, delt, L ,eqn);
[U3a, ~] = LF(a, CFLA, tspana, delt, L ,eqn);
[U4a, ~] = LW(a, CFLA, tspana, delt, L ,eqn);
[U5a, ~] = HighOrder(a, CFLA, tspana, delt, L ,eqn);

[U1b, ~] = FTBS(a, CFLA, tspanb, delt, L ,eqn);
[U2b, ~] = FTCS(a, CFLA, tspanb, delt, L ,eqn);
[U3b, ~] = LF(a, CFLA, tspanb, delt, L ,eqn);
[U4b, CFLA] = LW(a, CFLA, tspanb, delt, L ,eqn);
[U5b, ~] = HighOrder(a, CFLA, tspanb, delt, L ,eqn);
plotansCFLA = plotUsq2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
sgtitle('CFL = 0.50')
hold off

[U1a, ~] = FTBS(a, CFLB, tspana, delt, L ,eqn);
[U2a, ~] = FTCS(a, CFLB, tspana, delt, L ,eqn);
[U3a, ~] = LF(a, CFLB, tspana, delt, L ,eqn);
[U4a, ~] = LW(a, CFLB, tspana, delt, L ,eqn);
[U5a, ~] = HighOrder(a, CFLB, tspana, delt, L ,eqn);

[U1b, ~] = FTBS(a, CFLB, tspanb, delt, L ,eqn);
[U2b, ~] = FTCS(a, CFLB, tspanb, delt, L ,eqn);
[U3b, ~] = LF(a, CFLB, tspanb, delt, L ,eqn);
[U4b, CFLB] = LW(a, CFLB, tspanb, delt, L ,eqn);
[U5b, ~] = HighOrder(a, CFLB, tspanb, delt, L ,eqn);
plotansCFLB = plotUsq2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
sgtitle('CFL = 0.75')
hold off

[U1a, ~] = FTBS(a, CFLC, tspana, delt, L ,eqn);
[U2a, ~] = FTCS(a, CFLC, tspana, delt, L ,eqn);
[U3a, ~] = LF(a, CFLC, tspana, delt, L ,eqn);
[U4a, ~] = LW(a, CFLC, tspana, delt, L ,eqn);
[U5a, ~] = HighOrder(a, CFLC, tspana, delt, L ,eqn);

[U1b, ~] = FTBS(a, CFLC, tspanb, delt, L ,eqn);
[U2b, ~] = FTCS(a, CFLC, tspanb, delt, L ,eqn);
[U3b, ~] = LF(a, CFLC, tspanb, delt, L ,eqn);
[U4b, CFLC] = LW(a, CFLC, tspanb, delt, L ,eqn);
[U5b, ~] = HighOrder(a, CFLC, tspanb, delt, L ,eqn);
plotansCFLC = plotUsq2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
sgtitle('CFL = 0.98')
hold off








%% ACCURACY CALCULATION AND PLOT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section calls the accuracy.m function (in which already has loaded
% values) to calculate the L2 error for each algorithm
clear variables
close all
clc



% Calculating Accuracy for Square function:
[errors1,~,delt_values1] = accuracy(@FTBS, 1);
[errors2,~,delt_values2] = accuracy(@FTCS, 1);
[errors3,~,delt_values3] = accuracy(@LF, 1);
[errors4,~,delt_values4] = accuracy(@LW, 1);
[errors5,CFL_actual1,delt_values5] = accuracy(@HighOrder, 1);

% Plot Accuracy
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


% Calculating Accuracy for Sine-Squared function:
[errors1,~,delt_values1] = accuracy(@FTBS, 2);
[errors2,~,delt_values2] = accuracy(@FTCS, 2);
[errors3,~,delt_values3] = accuracy(@LF, 2);
[errors4,~,delt_values4] = accuracy(@LW, 2);
[errors5,CFL_actual2,delt_values5] = accuracy(@HighOrder, 2);

% Plot Accuracy
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





%% IDEAL MHD MODEL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section calculates the profiles of the state variables. There are
% two main sections, each are identical except for the tspan parameter.
% They calculate the variables and plot. tspan = 0.1 and 0.075
clear variables;
clc;
close all;

%%%%%%%%%%%%%%%%%% For tspan = 0.1 %%%%%%%%%%%%%%
% Defining initial conditions (Brio-Wu Shock Tube)
rho_l = 1;     u_l = 0; v_l = 0;  p_l = 1;    B_yl = 1;
rho_r = 0.125; u_r = 0; v_r = 0;  p_r = 0.1;  B_yr = -1;
B_x = 0.75;    gamma = 2;

e_l = 0.5*rho_l*(u_l^2+v_l^2) + p_l/(gamma-1) + 0.5*(B_x^2+B_yl^2);
e_r = 0.5*rho_r*(u_r^2+v_r^2) + p_r/(gamma-1) + 0.5*(B_x^2+B_yr^2);

Q_l = [rho_l; rho_l*u_l; rho_l*v_l; B_x; B_yl; e_l];
Q_r = [rho_r; rho_r*u_r; rho_r*v_r; B_x; B_yr; e_r];

% Parameter setup
L = 1;
J = 400;                 
x = linspace(0, L, J);
delx = x(2) - x(1);

% Initial condition split at L/2
U_q = zeros(6, J);
U_q(:, x < L/2) = repmat(Q_l, 1, sum(x < L/2));
U_q(:, x >= L/2) = repmat(Q_r, 1, sum(x >= L/2));

CFL = 0.4;

% Estimate max wavespeed (approx. fast magnetosonic speed)
a = sqrt(gamma * p_l / rho_l);
b_x = B_x / sqrt(rho_l);
Cfm = sqrt((a^2 + sqrt(a^4 + 4*b_x^2)) / 2);
delt = CFL * delx / Cfm;

tspan = 0.1;
N = floor(tspan / delt);
fprintf('Time step: %.5f, N steps: %d\n', delt, N);

% Time iteration loop
for n = 1:N
    U_f = zeros(6, J);

    % Compute flux at each point
    for j = 1:J
        rho = U_q(1,j);
        u   = U_q(2,j)/rho;
        v   = U_q(3,j)/rho;
        By  = U_q(5,j);
        e   = U_q(6,j);
        p   = (gamma - 1)*(e - 0.5*rho*(u^2 + v^2) - 0.5*(B_x^2 + By^2));
        
        U_f(:,j) = [ ...
            U_q(2,j);                                       % rho*u
            rho*u^2 + p + 0.5*(B_x^2 + By^2) - B_x^2;     % momentum-x
            rho*u*v - B_x*By;                             % momentum-y
            0;                                            % Bx (constant)
            u*By - v*B_x;                                 % By flux
            (e + p + 0.5*(B_x^2 + By^2))*u - B_x*(u*B_x + v*By) % energy
        ];
    end

    % Lax-Friedrichs scheme
    U_new = U_q;
    for j = 2:J-1
        U_new(:,j) = 0.5*(U_q(:,j+1) + U_q(:,j-1)) ...
                     - (delt/(2*delx)) * (U_f(:,j+1) - U_f(:,j-1));
    end

    % Reflective boundary conditions
    U_new(:,1) = U_q(:,2);
    U_new(:,end) = U_q(:,end-1);

    U_q = U_new;
end

% Recover primitive variables
rho = U_q(1,:);
u   = U_q(2,:)./U_q(1,:);
v   = U_q(3,:)./U_q(1,:);
By  = U_q(5,:);
e   = U_q(6,:);
p   = (gamma - 1) * (e - 0.5*rho.*(u.^2 + v.^2) - 0.5*(B_x^2 + By.^2));

% Plot results
figure(1);
subplot(2,3,1); plot(x, rho); title('Density \rho');
hold on
subplot(2,3,2); plot(x, u); title('Velocity u');
subplot(2,3,3); plot(x, v); title('Velocity v');
subplot(2,3,4); plot(x, p); title('Pressure p');
subplot(2,3,5); plot(x, By); title('B_y');
sgtitle('Brio & Wu MHD Solution (Lax-Friedrichs), tspan = 0.10');
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    grid(axesHandles(k), 'on');
    xlabel(axesHandles(k), 'x');
end
hold off










%%%%%%%%%%%%%% For tspan = 0.075; %%%%%%%%%%%%%
clear variables;


% Defining initial conditions (Brio-Wu Shock Tube)
rho_l = 1;     u_l = 0; v_l = 0;  p_l = 1;    B_yl = 1;
rho_r = 0.125; u_r = 0; v_r = 0;  p_r = 0.1;  B_yr = -1;
B_x = 0.75;    gamma = 2;

e_l = 0.5*rho_l*(u_l^2+v_l^2) + p_l/(gamma-1) + 0.5*(B_x^2+B_yl^2);
e_r = 0.5*rho_r*(u_r^2+v_r^2) + p_r/(gamma-1) + 0.5*(B_x^2+B_yr^2);

Q_l = [rho_l; rho_l*u_l; rho_l*v_l; B_x; B_yl; e_l];
Q_r = [rho_r; rho_r*u_r; rho_r*v_r; B_x; B_yr; e_r];

% Parameter setup
L = 1;
J = 400;                 
x = linspace(0, L, J);
delx = x(2) - x(1);

% Initial condition split at L/2
U_q = zeros(6, J);
U_q(:, x < L/2) = repmat(Q_l, 1, sum(x < L/2));
U_q(:, x >= L/2) = repmat(Q_r, 1, sum(x >= L/2));

CFL = 0.4;

% Estimate max wavespeed (approx. fast magnetosonic speed)
a = sqrt(gamma * p_l / rho_l);
b_x = B_x / sqrt(rho_l);
Cfm = sqrt((a^2 + sqrt(a^4 + 4*b_x^2)) / 2);
delt = CFL * delx / Cfm;

tspan = 0.075;
N = floor(tspan / delt);
fprintf('Time step: %.5f, N steps: %d\n', delt, N);

% Time iteration loop
for n = 1:N
    U_f = zeros(6, J);

    % Compute flux at each point
    for j = 1:J
        rho = U_q(1,j);
        u   = U_q(2,j)/rho;
        v   = U_q(3,j)/rho;
        By  = U_q(5,j);
        e   = U_q(6,j);
        p   = (gamma - 1)*(e - 0.5*rho*(u^2 + v^2) - 0.5*(B_x^2 + By^2));
        
        U_f(:,j) = [ ...
            U_q(2,j);                                       % rho*u
            rho*u^2 + p + 0.5*(B_x^2 + By^2) - B_x^2;     % momentum-x
            rho*u*v - B_x*By;                             % momentum-y
            0;                                            % Bx (constant)
            u*By - v*B_x;                                 % By flux
            (e + p + 0.5*(B_x^2 + By^2))*u - B_x*(u*B_x + v*By) % energy
        ];
    end

    % Lax-Friedrichs scheme
    U_new = U_q;
    for j = 2:J-1
        U_new(:,j) = 0.5*(U_q(:,j+1) + U_q(:,j-1)) ...
                     - (delt/(2*delx)) * (U_f(:,j+1) - U_f(:,j-1));
    end

    % Reflective boundary conditions
    U_new(:,1) = U_q(:,2);
    U_new(:,end) = U_q(:,end-1);

    U_q = U_new;
end

% Recover primitive variables
rho = U_q(1,:);
u   = U_q(2,:)./U_q(1,:);
v   = U_q(3,:)./U_q(1,:);
By  = U_q(5,:);
e   = U_q(6,:);
p   = (gamma - 1) * (e - 0.5*rho.*(u.^2 + v.^2) - 0.5*(B_x^2 + By.^2));

% Plot results
figure(1);
subplot(2,3,1); plot(x, rho); title('Density \rho');
hold on
subplot(2,3,2); plot(x, u); title('Velocity u');
subplot(2,3,3); plot(x, v); title('Velocity v');
subplot(2,3,4); plot(x, p); title('Pressure p');
subplot(2,3,5); plot(x, By); title('B_y');
sgtitle('Brio & Wu MHD Solution (Lax-Friedrichs), tspan = 0.075');
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    grid(axesHandles(k), 'on');
    xlabel(axesHandles(k), 'x');
end
hold off