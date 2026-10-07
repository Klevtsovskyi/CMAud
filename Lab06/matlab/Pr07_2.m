function Pr07_2
    clear all, close all, clc
    % вхiднi данi
    a = 0; b = 0.35;
    % параметри сітки і ітераційного процесу
    n = 13; eps = 1e-8; mi = 500;
    Y = 0.05.*ones(n,1);  % початкове наближення
    [X,Y] = frur_mk(a,b,1,n,@frurmk_k,@frurmk_f,Y,eps,mi,2);
    plot(X,Y,'k','LineWidth',2)
    title('Наближений розв''язок IP')
    xlabel('x'), ylabel('y'), grid on
end
function [y, yp] = frurmk_k(x, t, u)
    y = exp(x - u);  % опис ядра K(x,t,u)
    yp = -y;         % опис похідної dK(x,t,u)/du
end
function [y, yp] = frurmk_f(x, u)
    y = x - u;  % опис f(x,u)
    yp = -1;    % опис похідної df(x,u)/du
end