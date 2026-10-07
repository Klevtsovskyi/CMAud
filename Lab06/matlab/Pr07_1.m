% Pr07_1.m
clear all, close all, clc
% вхiднi данi
fr2mk_k = @(x,s)(1./(x+1));            % K(x,s)
fr2mk_f = @(x)(x.*x.*(1-x)-1./(x+1));  % f(x)
lamb = 12; a = 0; b = 1; n = 101;
X = linspace(a, b, n);
[x,YY] = fr2_mcol(lamb, X, fr2mk_k, fr2mk_f);
[X,Y,n] = fr2_mk(lamb, a, b, 1, n, fr2mk_k, fr2mk_f);
plot(X,Y,'m -',x,YY,'b :')
title('Наближений розв''язок IP y')
legend('м.КФ','м.колок','Location','Best')
xlabel('x'), ylabel('y'), grid on, pause, close all
u = @(x)(x.*x.*(1-x));    % u(x) - точний розв'язок
subplot(2,1,1); plot(X,u(X)-Y','r','LineWidth',2)
title('м.КФ: рiзниця розв''язків u-y')
xlabel('x'), ylabel('u-y'), grid on
subplot(2,1,2); plot(x,u(x)-YY','r','LineWidth',2)
title('м.Колокацій: рiзниця розв''язків u-y')
xlabel('x'), ylabel('u-y'), grid on