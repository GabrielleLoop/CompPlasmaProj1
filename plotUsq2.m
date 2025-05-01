function plotans = plotUsq2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
    % need to follow function with title and hold off
    % This function plots all algorithm results in relation to the Sine-Squared Algorithm

    [n1a,j1a] = size(U1a);
    xspan1a = linspace(0,L,j1a);
    eqn = @(x) (x >= pi/4 & x <= pi/2);
    realsol1a = eqn(xspan1a); %

    [n2a,j2a] = size(U2a);
    xspan2a = linspace(0,L,j2a);



    [n1b,j1b] = size(U1b);
    xspan1b = linspace(0,L,j1b);
    realsol1b = eqn(xspan1b); %

    [n2b,j1b] = size(U2b);
    xspan2b = linspace(0,L,j1b);
    [n3b,j3b] = size(U5b);
    xspan3b = linspace(0,L,j3b);



    Uplot1a = U1a(n1a,:);
    Uplot2a = U2a(n2a,:);
    Uplot3a = U3a(n2a,:);
    Uplot4a = U4a(n2a,:);
    Uplot5a = U5a(n1a,:);

    Uplot1b = U1b(n1b,:);
    Uplot2b = U2b(n2b,:);
    Uplot3b = U3b(n2b,:);
    Uplot4b = U4b(n2b,:);
    Uplot5b = U5b(n1b,:);



    plotans = figure;
    ax1 = subplot(2,1,1);
    pos1 = get(ax1,'Position');

    plot(xspan1a, Uplot1a)
    hold on
    plot(xspan2a, Uplot2a)
    plot(xspan2a, Uplot3a)
    plot(xspan2a, Uplot4a)
    plot(xspan1a, Uplot5a)
    plot(xspan1a,realsol1a, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('FTBS','FTCS','LF','LW','High Order', 'Real Solution')
   % legend('FTBS','LF','LW','High Order', 'Real Solution') % use if emmitting FTCS
    t1 = title('1 Period');
    t1.Position(1) = t1.Position(1) - 1.4;
    grid on
    hold off


    ax2 = subplot(2,1,2);
    pos2 = get(ax2, 'Position');
    
    plot(xspan1b, Uplot1b)
    hold on
    plot(xspan2b, Uplot2b)
    plot(xspan2b, Uplot3b)
    plot(xspan2b, Uplot4b)
    plot(xspan3b, Uplot5b)
    plot(xspan1b,realsol1b, 'LineStyle','--')
    xlabel('x')
    ylabel('U')
    legend('FTBS','FTCS','LF','LW', 'High Order','Real Solution')
   % legend('FTBS','LF','LW', 'High Order', 'Real Solution') % use if emmitting FTCS
    t2 = title('10 Periods');
    t2.Position(1) = t2.Position(1) - 1.4;
    grid on



    commonWidth = max(pos1(3), pos2(3));
    pos1(3) = commonWidth;
    pos2(3) = commonWidth;

    % pos2(2) = pos1(2) + commonWidth + 0.05;
    set(ax1, 'Position', pos1);
    set(ax2, 'Position', pos2);


end
