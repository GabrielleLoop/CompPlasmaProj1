% Forward Time Backward Space (FTBS), Forward Euler
function [U,CFL] = FTBS(a, CFL1, tspan, delt, L, eqn)
    % a = linear advection coefficient / wave speed
    % CFL1 = desired CFL number to analyze
    % tspan = total time to differentiate over
    % delt = differentiable time step
    % L = space domain
    % eqn = initial conditions of interest
    
    [CFL, J] = CFLdes(CFL1, a, delt, L);
    N = round(tspan/delt);
    x = linspace(0,L,J);
    
    % Using functions to Solve for Initial Conditions:
    U = zeros(N,J);
    U(1,:) = eqn(x);
    u_new = zeros(1,J);
    
    
    % Populating the rest of U:
    for n = 1:N-1
        for j = 2:J
            u_new(:,j) = U(n,j)- CFL * (U(n,j)-U(n,j-1));
        end
        u_new(:,1) = u_new(:,J);
        U(n+1,:) = u_new;
        u_new = zeros(1,J);
    end


end
