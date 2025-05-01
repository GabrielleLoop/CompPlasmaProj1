clear variables;
clc;
close all;

% === INITIAL CONDITIONS (Brio-Wu Shock Tube) ===
rho_l = 1;     u_l = 0; v_l = 0;  p_l = 1;    B_yl = 1;
rho_r = 0.125; u_r = 0; v_r = 0;  p_r = 0.1;  B_yr = -1;
B_x = 0.75;    gamma = 2;

e_l = 0.5*rho_l*(u_l^2+v_l^2) + p_l/(gamma-1) + 0.5*(B_x^2+B_yl^2);
e_r = 0.5*rho_r*(u_r^2+v_r^2) + p_r/(gamma-1) + 0.5*(B_x^2+B_yr^2);

Q_l = [rho_l; rho_l*u_l; rho_l*v_l; B_x; B_yl; e_l];
Q_r = [rho_r; rho_r*u_r; rho_r*v_r; B_x; B_yr; e_r];

% === GRID SETUP ===
L = 1;
J = 400;                  % Number of spatial points
x = linspace(0, L, J);
delx = x(2) - x(1);

% Initial condition split at L/2
U_q = zeros(6, J);
U_q(:, x < L/2) = repmat(Q_l, 1, sum(x < L/2));
U_q(:, x >= L/2) = repmat(Q_r, 1, sum(x >= L/2));

% === TIME STEP SETUP ===
CFL = 0.4;

% Estimate max wavespeed (approx. fast magnetosonic speed)
a = sqrt(gamma * p_l / rho_l);
b_x = B_x / sqrt(rho_l);
Cfm = sqrt((a^2 + sqrt(a^4 + 4*b_x^2)) / 2);
delt = CFL * delx / Cfm;

tspan = 0.1;
N = floor(tspan / delt);
fprintf('Time step: %.5f, N steps: %d\n', delt, N);

% === MAIN TIME INTEGRATION LOOP ===
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

    % Lax-Friedrichs scheme (simple Riemann approximation)
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

% === RECOVER PRIMITIVE VARIABLES FOR PLOTTING ===
rho = U_q(1,:);
u   = U_q(2,:)./U_q(1,:);
v   = U_q(3,:)./U_q(1,:);
By  = U_q(5,:);
e   = U_q(6,:);
p   = (gamma - 1) * (e - 0.5*rho.*(u.^2 + v.^2) - 0.5*(B_x^2 + By.^2));

% === PLOT RESULTS ===
figure(1);
subplot(2,3,1); plot(x, rho); title('Density \rho');
hold on
subplot(2,3,2); plot(x, u); title('Velocity u');
subplot(2,3,3); plot(x, v); title('Velocity v');
subplot(2,3,4); plot(x, p); title('Pressure p');
subplot(2,3,5); plot(x, By); title('B_y');
% subplot(2,3,6); plot(x, e); title('Total Energy');
sgtitle('Brio & Wu MHD Solution (Lax-Friedrichs), tspan = 0.10');
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    grid(axesHandles(k), 'on');
    xlabel(axesHandles(k), 'x');
end
hold off










%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;
%clc;
%close all;

% === INITIAL CONDITIONS (Brio-Wu Shock Tube) ===
rho_l = 1;     u_l = 0; v_l = 0;  p_l = 1;    B_yl = 1;
rho_r = 0.125; u_r = 0; v_r = 0;  p_r = 0.1;  B_yr = -1;
B_x = 0.75;    gamma = 2;

e_l = 0.5*rho_l*(u_l^2+v_l^2) + p_l/(gamma-1) + 0.5*(B_x^2+B_yl^2);
e_r = 0.5*rho_r*(u_r^2+v_r^2) + p_r/(gamma-1) + 0.5*(B_x^2+B_yr^2);

Q_l = [rho_l; rho_l*u_l; rho_l*v_l; B_x; B_yl; e_l];
Q_r = [rho_r; rho_r*u_r; rho_r*v_r; B_x; B_yr; e_r];

% === GRID SETUP ===
L = 1;
J = 400;                  % Number of spatial points
x = linspace(0, L, J);
delx = x(2) - x(1);

% Initial condition split at L/2
U_q = zeros(6, J);
U_q(:, x < L/2) = repmat(Q_l, 1, sum(x < L/2));
U_q(:, x >= L/2) = repmat(Q_r, 1, sum(x >= L/2));

% === TIME STEP SETUP ===
CFL = 0.4;

% Estimate max wavespeed (approx. fast magnetosonic speed)
a = sqrt(gamma * p_l / rho_l);
b_x = B_x / sqrt(rho_l);
Cfm = sqrt((a^2 + sqrt(a^4 + 4*b_x^2)) / 2);
delt = CFL * delx / Cfm;

tspan = 0.075;
N = floor(tspan / delt);
fprintf('Time step: %.5f, Nsteps: %d\n', delt, N);

% === MAIN TIME INTEGRATION LOOP ===
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

    % Lax-Friedrichs scheme (simple Riemann approximation)
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

% === RECOVER PRIMITIVE VARIABLES FOR PLOTTING ===
rho = U_q(1,:);
u   = U_q(2,:)./U_q(1,:);
v   = U_q(3,:)./U_q(1,:);
By  = U_q(5,:);
e   = U_q(6,:);
p   = (gamma - 1) * (e - 0.5*rho.*(u.^2 + v.^2) - 0.5*(B_x^2 + By.^2));

% === PLOT RESULTS ===
figure(2);
subplot(2,3,1); plot(x, rho); title('Density \rho');
hold on
subplot(2,3,2); plot(x, u); title('Velocity u');
subplot(2,3,3); plot(x, v); title('Velocity v');
subplot(2,3,4); plot(x, p); title('Pressure p');
subplot(2,3,5); plot(x, By); title('B_y');
% subplot(2,3,6); plot(x, e); title('Total Energy');
sgtitle('Brio & Wu MHD Solution (Lax-Friedrichs), tspan = 0.075');
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    grid(axesHandles(k), 'on');
    xlabel(axesHandles(k), 'x');
end
hold off

