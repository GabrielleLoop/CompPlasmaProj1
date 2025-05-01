% Lax-Wendroff
function [U,CFL] = LW(a, CFL1, tspan, delt, L, eqn)
    % a = linear advection coefficient
    % tspan = total time to differentiate over
    % L = space domain
    % N = number of time points or iterations
    % J = number of space points or iterations
    % U_j1, U_j2, U_j3 = initial conditions U(1,j-1), U(1,j), U(1,j+1)
    
    [CFL, J] = CFLdes(CFL1, a, delt, L);
    N = tspan/delt;
    x = linspace(0,L,J+1);
    
    % Using functions to Solve for Initial Conditions:
    U = zeros(N,J+1);
    U(1,:) = eqn(x);
    u_new = zeros(1,J+1);
    
    
    % Populating U:
    for n=1:N
        for j = 2:J
        % Step 1: Lax-Friedrichs:
        u_new_half_for = (U(n,j+1)+U(n,j))/2 - CFL/2*(U(n,j+1)-U(n,j)); % (:,j)
        u_new_half_back = (U(n,j)+U(n,j-1))/2 - CFL/2*(U(n,j)-U(n,j-1)); % (:.j-1)
        % Step 2: Leap Frog:
        u_new(:,j) = U(n,j) - CFL*(u_new_half_for - u_new_half_back);
        end
        u_new(:,1) = u_new(:,J);
        U(n+1,:) = u_new;
        u_new = zeros(1,J+1);
    end

end

