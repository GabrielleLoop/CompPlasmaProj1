function [U, CFL] = HighOrder(a, CFL1, tspan, delt, L, eqn)
    % a     = advection speed
    % CFL1  = target CFL number
    % tspan = final time
    % delt  = time step
    % L     = domain length
    % eqn   = initial condition function handle

    [CFL, J] = CFLdes(CFL1, a, delt, L);  % Get actual CFL and J
    N = round(tspan / delt);             % Number of time steps
    x = linspace(0, L, J);               % Spatial grid
    delx = x(2) - x(1);

    % Initial condition
    U = zeros(N, J);
    U(1,:) = eqn(x);
    dU1 = zeros(1,J);
    dU2 = zeros(1,J);

    for n = 1:N
        % Compute spatial derivative using 4th order centered difference
        for j = 1:J
            jm2 = mod(j-3, J) + 1;
            jm1 = mod(j-2, J) + 1;
            jp1 = mod(j, J) + 1;
            jp2 = mod(j+1, J) + 1;
    
            dU1(j) = a/(12*delx) * (-U(jp2) + 8*U(jp1) - 8*U(jm1) + U(jm2));
        end
        % Stage 1: Euler step
        U_1a = U(n,:) - delt * dU1;

        % Stage 2: compute spatial derivative at U_star
        for j = 1:J
            jm2 = mod(j - 3, J) + 1;
            jm1 = mod(j - 2, J) + 1;
            jp1 = mod(j    , J) + 1;
            jp2 = mod(j + 1, J) + 1;
    
            dU2(j) = a/(12*delx) * (-U_1a(jp2) + 8*U_1a(jp1) - 8*U_1a(jm1) + U_1a(jm2));
        end

        % RK2 step
        if n~=N
        U(n+1,:) = U(n,:) - 0.5 * delt * (dU1 + dU2);
        U(:,1) = U(:,J);
        end
    end
end

% =========================
% 4th-order centered space derivative
% Periodic BCs assumed
%function dU = spatial_derivative_4th_order(U, dx, a)
%    J = length(U);
%    dU = zeros(1, J);
%
%    for j = 1:J
%        jm2 = mod(j - 3, J) + 1;
%        jm1 = mod(j - 2, J) + 1;
%        jp1 = mod(j    , J) + 1;
%        jp2 = mod(j + 1, J) + 1;
%
%        dU(j) = a/(12*dx) * (-U(jp2) + 8*U(jp1) - 8*U(jm1) + U(jm2));
%    end
%end