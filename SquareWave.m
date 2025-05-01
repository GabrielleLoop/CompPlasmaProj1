%%%%%%%%%%%%%%%%% Square Wave %%%%%%%%%%%%%%%%%%%%%%%%%
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











%%
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


