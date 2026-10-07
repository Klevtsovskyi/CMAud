function Pr07_5
    clear all, close all, clc
    lamb = 0.25; a = 0; b = 1;   % вхiднi данi
    n = 13; eps = 1e-7; mi = 5000;
    [X,Y,n] = fr2_mpn(lamb,a,b,1,n,eps,mi,@fr2mpn_k,@fr2mpn_f);
    plot(X,Y,'m','LineWidth',2)
    title('Наближений розв''язок IP y')
    xlabel('x'), ylabel('y'), grid on, pause, close all
    u = @(x)(x.*x);    % u(x) - точний розв'язок
    Z = u(X)-Y';       % Рiзниця мiж точним i набл.розв'язком
    plot(X,Z,'r','LineWidth',2)
    title('Рiзниця розв''язків u-y')
    xlabel('x'), ylabel('u-y'), grid on
end
function yp = fr2mpn_k(x,s)
    yp = x.*s+1;            %   K(x,s)
end
function yp = fr2mpn_f(x)
    yp = (x-1/16).*x-1/12;  %   f(x)
end