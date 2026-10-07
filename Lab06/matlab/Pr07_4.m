function Pr07_4
    clear all, close all, clc
    %  а) модельна задача
    lamb = -1; a = 0; b = 2;        % вхiднi данi
    n = 51;
    [X,Y,n] = vo2_mk(lamb, a, b, n, @vo2mk_k, @vo2mk_f);
    plot(X,Y,'m','LineWidth',2) 
    title('Наближений розв''язок IP y')
    xlabel('x'), ylabel('y'), grid on, pause, close all
    u = @(x)(1./sqrt((1+x.*x).^3));    % u(x) - точний розв'язок
    Z = u(X)-Y'; % Рiзниця мiж точним i набл.розв'язком
    plot(X,Z,'r'), title('Рiзниця розв''язків u-y')
    xlabel('x'), ylabel('u-y'), grid on, pause, close all

    %  б)
    lamb = 1; a = 0; b = 2;         % вхiднi данi
    [X,Y,n] = vo2_mk(lamb, a, b, n, @vo2_k, @vo2_f);
    plot(X,Y,'k','LineWidth',2), grid on
    title('Наближений розв''язок IP y')
    xlabel('x'), ylabel('y')
end
function yp = vo2mk_k(x,s)
    yp = s./(1+x.*x);  %  K(x,s)
end
function yp = vo2mk_f(x)
    yp = 1./(1+x.*x);  %  f(x)
end
function yp = vo2_k(x,s)
    yp = x.*s-1;       %  K(x,s)
end
function yp = vo2_f(x)
    yp = x;            %  f(x)
end