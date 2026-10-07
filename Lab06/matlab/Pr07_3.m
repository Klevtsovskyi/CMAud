function Pr07_3
    clear all, close all, clc
    % a) модельна задача
    global m
    lamb = 1/3; a = 0; b = 1;   % вхiднi данi
    n = 21;   % кількість кв.вузлів м.Сімпсона
    m = 2;    
    [X,Y,n] = fr2_mzjd(lamb,a,b,1,n,@fr2mzj_k,@fr2mzj_f);
    plot(X,Y,'m','LineWidth',2) 
    title('Наближений розв''язок IP y')
    xlabel('x'), ylabel('y'), grid on, pause, close all
    u = @(x)(x);    % u(x) - точний розв'язок
    Z = u(X)-Y';    % Рiзниця мiж точним i набл.розв'язком
    plot(X,Z,'r','LineWidth',2)
    title('Рiзниця розв''язків u-y')
    xlabel('x'), ylabel('u-y'), grid on, pause, close all
    % б)
    lamb = 1; a = 0; b = 2;     % вхiднi данi
    cl = 'mbkr'; M = [4:4:16];
    for j = 1:4;
        m = M(j); L{j} = sprintf('m=%u',m);
        [X,Y,n] = fr2_mzjd(lamb, a, b, 1, n, @mzj_k, @mzj_f);
        plot(X,Y,cl(j)), hold on
    end
    hold off, title('Наближений розв''язок IP')
    grid on, xlabel('x'), ylabel('y')
    legend(L{:},'Location','Best')
end
function [A, B] = fr2mzj_k(x)
    % ядро K(x,s)=SUM(k=1,m,Ak(x)*Bk(s))
    global m
    t = x.^2;
    A = [x; t];
    B = [t; x];
end
function yp = fr2mzj_f(x)
    yp = (11/12-1/9.*x).*x;  %  f(x)
end
function [A,B] = mzj_k(x)
    % опис ядра  K(x,s)=SUM(k=1,m,Ak(x)*Bk(s))
    global m
    f = 1; z = 1;
    A = zeros(m,1); B = A;
    for i = 1:m
       z = z.*x; f = f.*i;
       A(i) = z; B(i) = z./f;
    end
end
function yp = mzj_f(x)
    yp = x;  % f(x)
end