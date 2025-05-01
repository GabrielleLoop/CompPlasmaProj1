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



