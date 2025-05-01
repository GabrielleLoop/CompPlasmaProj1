clear var
close all
clc




%% Square Wave
% u(x) = if pi/4
clc
clear variables
close all


% Defining equation:
eqn = @(x) (x >= pi/4 & x <= pi/2);

% Defining solving parameters:
L = pi;
J = 101;
% N = 200;
a = 1;
CFL = 0.97;

delx = L/(J+1);
delt = CFL*delx/a;
tspan = 1*pi; 
N = round(tspan/delt);

CFLa = 0.35;
delta = CFLa*delx/a;
Na = round(tspan/delta);

CFLb = 0.97;
deltb = CFLb*delx/a;
Nb = round(tspan/deltb);

% Testing different algorithms:
% Forward Time Backward Space:
[Ua,CFLa,~] = FTBS(a, CFLa, Na, L, J, eqn);
[Ub,CFLb,~] = FTBS(a, CFLb, Nb, L, J, eqn);
plotUsq(L, Ua(N,:));
title('FTBS Sine-Squared Wave, CFL = 0.35')
hold off
plotUsq(L, Ub(N,:))
title('FTBS Sine-Squared Wave, CFL = 0.97')
hold off

% Forward Time Center Space:
U2a = FTCS(a, CFLa, Na, L, J, eqn);
U2b = FTCS(a, CFLb, Nb, L, J, eqn);
plotUsq(L, U2a(N,:));
title('FTCS Sine-Squared Wave, CFL = 0.35')
hold off
plotUsq(L, U2b(N,:))
title('FTCS Sine-Squared Wave, CFL = 0.97')
hold off

% Lax-Friedrichs:
U3a = LF(a, CFLa, Na, L, J, eqn);
U3b = LF(a, CFLb, Nb, L, J, eqn);
plotUsq(L, U3a(N,:));
title('LF Sine-Squared Wave, CFL = 0.35')
hold off
plotUsq(L, U3b(N,:))
title('LF Sine-Squared Wave, CFL = 0.97')
hold off

% Lax-Wendroff
U4a = LW(a, CFLa, Na, L, J, eqn);
U4b = LW(a, CFLb, Nb, L, J, eqn);
plotUsq(L, U4a(N,:));
title('Lax-Wendroff Sine-Squared Wave, CFL = 0.35')
hold off
plotUsq(L, U4b(N,:))
title('Lax-Wendroff Sine-Squared Wave, CFL = 0.97')
hold off


% High-Order 4-2



















%% Sine-Squared Wave
% u(x) = sin^2(kx), k=1
clc
clear variables
close all


% Defining equation:
k = 1;
eqn = @(x) (sin(k.*x)).^2;


% Defining solving parameters:
L = pi;
J = 101;
% N = 200;
a = 1;
CFL = 0.97;

delx = L/(J+1);
delt = CFL*delx/a;
tspan = 1*pi; 
N = round(tspan/delt);

CFLa = 0.35;
delta = CFLa*delx/a;
Na = round(tspan/delta);

CFLb = 0.97;
deltb = CFLb*delx/a;
Nb = round(tspan/deltb);

% Testing different algorithms:
% Forward Time Backward Space:
[Ua,CFLa,~] = FTBS(a, CFLa, Na, L, J, eqn);
[Ub,CFLb,~] = FTBS(a, CFLb, Nb, L, J, eqn);
plotUsin(L, Ua(N,:));
title('FTBS Sine-Squared Wave, CFL = 0.35')
hold off
plotUsin(L, Ub(N,:))
title('FTBS Sine-Squared Wave, CFL = 0.97')
hold off

% Forward Time Center Space:
U2a = FTCS(a, CFLa, Na, L, J, eqn);
U2b = FTCS(a, CFLb, Nb, L, J, eqn);
plotUsin(L, U2a(N,:));
title('FTCS Sine-Squared Wave, CFL = 0.35')
hold off
plotUsin(L, U2b(N,:))
title('FTCS Sine-Squared Wave, CFL = 0.97')
hold off

% Lax-Friedrichs:
U3a = LF(a, CFLa, Na, L, J, eqn);
U3b = LF(a, CFLb, Nb, L, J, eqn);
plotUsin(L, U3a(N,:));
title('LF Sine-Squared Wave, CFL = 0.35')
hold off
plotUsin(L, U3b(N,:))
title('LF Sine-Squared Wave, CFL = 0.97')
hold off

% Lax-Wendroff
U4a = LW(a, CFLa, Na, L, J, eqn);
U4b = LW(a, CFLb, Nb, L, J, eqn);
plotUsin(L, U4a(N,:));
title('Lax-Wendroff Sine-Squared Wave, CFL = 0.35')
hold off
plotUsin(L, U4b(N,:))
title('Lax-Wendroff Sine-Squared Wave, CFL = 0.97')
hold off


% High-Order 4-2






%%
% Defining solving parameters:
L = pi;
J = 101;
% N = 200;
a = 1;
CFL = 0.97;

delx = L/(J+1);
delt = CFL*delx/a;
tspan = pi; 
N = round(tspan/delt);


% Testing different algorithms:
% Forward Time Backward Space:
[U,CFL,x] = FTBS(a, CFL, N, L, J, eqn);
answerFTBS2 = plotUsin(L, U(N,:));
title('FTBS Sine-Squared Wave')
hold off

% Forward Time Center Space:
U2 = FTCS(a, CFL, N, L, J, eqn);
answerFTCS = plotUsin(L, U2(N,:));
title('FTCS Sine-Squared Wave')
hold off

% Lax-Friedrichs:
U3 = LF(a, CFL, N, L, J, eqn);
answerLF = plotUsin(L, U3(N,:));
title('Lax-Friedrichs Sine-Squared Wave')
hold off

% Lax-Wendroff
U4 = LW(a, CFL, N, L, J, eqn);
answerLW2 = plotUsin(L, U4(N,:));
title('Lax-Wendroff Sine-Squared Wave')
hold off


% High-Order 4-2









%% IDEAL MHD: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc

rho_l = 1;
u_l = 0;
v_l = 0;
p_l = 1;
B_yl = 1;
rho_r = 0.125;
u_r = 0;
v_r = 0;
p_r = 0.1;
B_yr = -1;
B_x = 0.75; % ???
g = 2;

e_l = 1/2*rho_l*(u_l^2+v_l^2) + p_l/(g-1) + 1/2*(B_x^2+B_yl^2);
e_r = 1/2*rho_r*(u_r^2+v_r^2) + p_r/(g-1) + 1/2*(B_x^2+B_yr^2);

% Defining equation:
Q_l = [rho_l; rho_l*u_l; rho_l*v_l; B_x; B_yl; e_l];
Q_r = [rho_r; rho_r*u_r; rho_r*v_r; B_x; B_yr; e_r];
eqnQ = @(x) (x > 10000) .* Q_r + (x < 10000) .* Q_l;

% initial data consists of two constant states U_1, U_r
% Q = U = [rho, rho*u, rho*v, B_x, B_y, e]^T
% e = 1/2*rho*u^2 + p/(g-1) + 1/2*B^2


F_l = [rho_l*u_l;
    rho_l*u_l^2 + p_l + 1/2*(B_x^2+B_yl^2) - B_x^2;
    rho_l*u_l*v_l - B_x*B_yl;
    0;
    u_l*B_yl - v_l*B_x;
    u_l*(e_l + p_l + 1/2*(B_x^2+B_yl^2)) - B_x*(u_l*B_x+v_l*B_yl)];
F_r = [rho_r*u_r;
    rho_r*u_r^2 + p_r + 1/2*(B_x^2+B_yr^2) - B_x^2;
    rho_r*u_r*v_r - B_x*B_yr;
    0;
    u_r*B_yr - v_r*B_x;
    u_r*(e_r + p_r + 1/2*(B_x^2+B_yr^2)) - B_x*(u_r*B_x+v_r*B_yr)];

eqnF = @(x) (x > 0.5) .* F_r + (x < 0.5) .* F_l;




a_star = sqrt((g*p_l + (B_x+B_yl))/(rho_l));
a = sqrt(g*p_l/rho_l);
b_x = B_x/sqrt(rho_l);
Cfm = sqrt((a_star^2+sqrt(a_star^4-4*a^2*b_x^2))/2);
CFL = 0.5;




% time steps = 10,000; grid pts = 20,000
% Defining solving parameters:
L = 1; % ?
J = 2000;
% tspan = 0.1;
N = 1000;




x = linspace(0,L,J+1);
delx = L/(J+1);
delt = CFL*delx/Cfm;
tspan = N*delt;




% Using functions to Solve for Initial Conditions:
U_q = zeros(6,J+1);
U_q(1:6,:) = eqnQ(x);
U_f = zeros(6,J+1);
U_f(1:6,:) = eqnF(x);
u_new = zeros(6,J+1);

% Populating U:
for n=1:N
    for j = 2:J
    u_new(:,j) = (U_q(:,j+1)+U_q(:,j-1))./2 - CFL/2.*(U_f(:,j+1)-U_f(:,j-1));
    end
    % u_new(:,j-1) = u_new(:,j);
    U(n*6-5:n*6,:) = u_new;
    u_new = zeros(6,J+1);
end

rho_final = U(1,:);
rhou = U(2,:);
u_final = rho_final./rhou;

rhov = U(3,:);
v_final = rho_final./rhov;

B_x_final = U(4,:);
B_y_final = U(5,:);

e_final = U(6:201:6,:);


figure(1)
subplot(5,1,1);
plot(x,u_final);
grid on

subplot(5,1,2);
plot(x,v_final);
grid on

subplot(5,1,3);
plot(x,B_x_final);
grid on

subplot(5,1,4);
plot(x,B_y_final);
grid on

subplot(5,1,5);
plot(x,e_final);
grid on




%%
% Using functions to Solve for Initial Conditions:
U_q = zeros(N*6,J+1);
U_q(1:6,:) = eqnQ(x);
U_f = zeros(N*6,J+1);
U_f(1:6,:) = eqnF(x);
u_new = zeros(6,J+1);

% Populating U:
for n=1:N
    for j = 2:J
    u_new(:,j) = (U_q(n,j+1)+U_q(n,j-1))./2 - CFL/2.*(U_f(n,j+1)-U_f(n,j-1));
    end
    % u_new(:,j-1) = u_new(:,j);
    U(n*6-5:n*6,:) = u_new;
    u_new = zeros(6,J+1);
end

rho_final = U(1:201-5:6,:);
rhou = U(2:201-4:6,:);
u_final = rho_final./rhou;

rhov = U(3:201-3:6,:);
v_final = rho_final./rhov;

B_x_final = U(4:201-2:6,:);
B_y_final = U(5:201-1:6,:);

e_final = U(6:201:6,:);




figure(1)
subplot(5,1,1);
plot(x,u_final);
grid on

subplot(5,1,2);
plot(x,v_final);
grid on

subplot(5,1,3);
plot(x,B_x_final);
grid on

subplot(5,1,4);
plot(x,B_y_final);
grid on

subplot(5,1,5);
plot(x,e_final);
grid on













%% %%%%%%%%%%%%%%%%% FUNCTIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting:
function plotans = plotUsin(L,U)
    % need to follow function with title and hold off

    [~,j] = size(U);
    xspan = linspace(0,L,j);
    k = 1;
    realsol = (sin(k.*xspan)).^2;

    Uplot1 = U;
    plotans = figure;
    plot(xspan, Uplot1)
    hold on
    plot(xspan,realsol, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('Algorithm', 'Real Solution')
    grid on


end

function plotans = plotUsin2(L,U1,U2)
    % need to follow function with title and hold off

    [~,j] = size(U1);
    xspan = linspace(0,L,j);
    k = 1;
    realsol = (sin(k.*xspan)).^2;

    Uplot1 = U1;
    Uplot2 = U2;
    plotans = figure;
    subplot(2,1,1);
    plot(xspan, Uplot1)
    hold on
    plot(xspan,realsol, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('Algorithm', 'Real Solution')
    grid on
    subplot(2,1,2);    
    plot(xspan, Uplot2)
    hold on
    plot(xspan,realsol, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('Algorithm', 'Real Solution')
    grid on


end


function plotans = plotUsq(L,U)
    % need to follow function with title and hold off

    [~,j] = size(U);
    xspan = linspace(0,L,j);
    eqn = @(x) (x >= pi/4 & x <= pi/2);
    realsol = eqn(xspan); % 
    Uplot1 = U;
    plotans = figure;
    plot(xspan, Uplot1)
    hold on
    plot(xspan,realsol, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('Algorithm', 'Real Solution')
    grid on


end


function plotans = plotUsq2(L,U1,U2)
    % need to follow function with title and hold off

    [~,j] = size(U1);
    xspan = linspace(0,L,j);
    eqn = @(x) (x >= pi/4 & x <= pi/2);
    realsol = eqn(xspan); %

    Uplot1 = U1;
    Uplot2 = U2;
    plotans = figure;
    plot(xspan, Uplot1)
    hold on
    plot(xspan,realsol, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('Algorithm', 'Real Solution')
    grid on
    hold off

    subplot(2,1,2);    
    plot(xspan, Uplot2)
    hold on
    plot(xspan,realsol, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('Algorithm', 'Real Solution')
    grid on
    hold off


end
































%% High-Order 4-2 (fourth-order space, second-order time)
