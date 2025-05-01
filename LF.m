% Lax-Friedrichs (LF)
function [U,CFL] = LF(a, CFL1, tspan, delt, L, eqn)
    % a = linear advection coefficient / wave speed
    % CFL1 = desired CFL number to analyze
    % tspan = total time to differentiate over
    % delt = differentiable time step
    % L = space domain
    % eqn = initial conditions of interest
    
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
        u_new(:,j) = (U(n,j+1)+U(n,j-1))/2 - CFL/2*(U(n,j+1)-U(n,j-1));
        end
        u_new(:,1) = u_new(:,J);
        U(n+1,:) = u_new;
        u_new = zeros(1,J+1);
    end

end
