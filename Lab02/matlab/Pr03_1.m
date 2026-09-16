function Pr03_1
% Наближене знаходження кореня f(x)=0

    close all, clear all, clc, warning off

    x1 = linspace(-2,4,161); y1 = psi(x1);
    x2 = linspace(-2,0,81); y2 = phi(x2);
    x3 = linspace(2,4,41); x3(1) = x3(1)+1e-2; y3 = phi(x3);

    plot(x1,y1,'b--','LineWidth',2); hold on
    plot(x2,y2,'k',x3,y3,'k','LineWidth',2); hold off; grid on
    legend('{\psi}','{\phi}','Location','Best');
    D = input('Діапазон кореня=');
    z = input('Наближення до кореня=');
    close all; clear x1 y1

    ep = 1e-8; mki = 100;  % ітераційні параметри

    global tau
    X = linspace(D(1),D(2),51); Y = df(X);
    tau = 1./(max(abs(Y))+ep);
    if df(z)>0
        tau = -tau;
    end
    S = ds(X); q = max(abs(S));
    plot(x2,f(x2),'r',x3,f(x3),'r','LineWidth',2);
    grid on; xlabel('x'); ylabel('f(x)');
    title('f(x)'); pause; close all
    disp(strcat('q=',num2str(q)));
    plot(X,S,'m-.','LineWidth',2);
    grid on; xlabel('x'); ylabel('s''(x)');
    title('s''(x)'); pause; close all
    clear x2 x3 y2 y3 X Y S

    fprintf('%-17s %21s %21s %5s %3s\n',...
            'Ітераційний метод','    X-наближен.корінь',...
            '    Значення f(x) в Х','КілІт','err');
    fmt = '%-17s %21.12e %21.12e %5u %3u\n';

    [x, k, err] = m_1dihot(@f, D, ep, mki);
    fprintf(fmt,'діхотомії        ',x,f(x),k,err);

    x = z; [x, k, err] = m_1simpit(@s, x, ep, mki);
    fprintf(fmt,'простої ітерації ',x,f(x),k,err);

    x = z; [x, k, err] = m_1estef(@s, x, ep, mki);
    fprintf(fmt,'Ейтк.-Стеффенсена',x,f(x),k,err);

    x = z; [x, k, err] = m_1njut(@f, x, ep, mki, @df, 1);
    fprintf(fmt,'Ньютона          ',x,f(x),k,err);

    x = z; x0 = x-0.05; [x, k, err] = m_1secant(@f,x,ep,mki,x0);
    fprintf(fmt,'січних           ',x,f(x),k,err);

    x = z; x0 = x-0.05; [x, k, err] = m_1kurch(@f,x,ep,mki,x0);
    fprintf(fmt,'Курчатова        ',x,f(x),k,err);
end

function y = phi(x)
% опис phi(x)
    y = sqrt(x.^3./(x-2));
end

function y = psi(x)
% опис psi(x)
    y = 13.1-1.5672.*(x-1.2).^2;
end

function y = f(x)
% опис f(x)
    y = phi(x)-psi(x);
end

function y = df(x)
% опис f'(x)
    y = 3.1344.*(x-1.2)+(x-3).*sqrt(x./(x-2).^3);
end

function y = s(x)
% опис s(x)
    global tau
    y = x+tau.*f(x);
end

function y = ds(x)
% опис s'(x)
    global tau
    y = 1.0+tau.*df(x);
end