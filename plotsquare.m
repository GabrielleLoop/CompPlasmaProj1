function [plotted1, plotted2, plotted3] = plotsquare(a, CFLA, tspana, tspanb, delt, L ,eqn)
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
[plotted1,plotted2,plotted3] = plotUsq2(L,U1a,U2a,U3a,U4a,U5a,U1b,U2b,U3b,U4b,U5b)
end