function [CFL_real, J_real] = CFLdes(CFL_des, a, delt, L)
% This function uses a desired CFL number and computes J and the CFL again
% so that J is an integer:
    delx_des = a*delt/CFL_des;
    J_des = L/delx_des - 1;
    J_real = round(J_des);
    CFL_real = a*delt/(L/J_real);
end