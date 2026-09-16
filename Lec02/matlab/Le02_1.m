function ns_main
% контрольний розрахунок для системи нелінійних рівнянь
%   F(X) = 0;  де Xточн = 1.

    close all, clear all, clc, warning off
    global fun

    N = 200; Xe = ones(N,1); X0 = 0.5.*Xe;

    e = 1e-10; mki = 200;  % параметри ітераційного процесу
    pr = 0;  % 1 - друк розв'язку і значення функції (0 - немає)

    disp(strcat('Точність ітераційн.процесу=', num2str(e)));
    disp(strcat('Обмеження кількості ітерац=', int2str(mki)));
    disp(strcat('Кількість невідомих=', int2str(N)));
    fprintf('%-18s %-8s %-8s %-20s\n','Метод  розв''язання',...
        'Кільк.іт','Час(сек)','Похибка');

    X = X0; t1 = clock;
    [X, it, err] = m_ns_njut(@FdF, X, e, mki); t2 = clock;
    print('Ньютона',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; t1 = clock;
    [X, it, err] = m_ns_njutg(@FdF, X, e, mki); t2 = clock;
    print('Ньютон+однов.пошук',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; t1 = clock;
    [X, it, err] = m_ns_broy(@FdF, X, e, mki); t2 = clock;
    print('Бройдена',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; t1 = clock;
    [X, it, err] = m_ns_pirs(@FdF, X, e, mki); t2 = clock;
    print('Пірсона',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; t1 = clock;
    [X, it, err] = m_ns_broym(@FdF, X, e, mki); t2 = clock;
    print('модифік. Бройдена',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; t1 = clock;
    [X, it, err] = m_ns_sym1r(@FdF, X, e, mki); t2 = clock;
    print('симетр.1-го рангу',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; X1 = X0 - 1e-3; t1 = clock;
    [X, it, err] = m_ns_kurch(@FdF, @Dd, X, e, mki, X1); t2 = clock;
    print('Курчатова',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

    X = X0; t1 = clock;
    [X, it, err] = m_ns_grad(@FdF, X, e, mki); t2 = clock;
    print('градієнтний',it,etime(t2,t1),norm(X-Xe,inf),X,pr);

end

function [varargout] = FdF(X)
% підфункція для обчислення f = {F(X)} при fun=0,
% а при fun=1 - пари {f, J}, де J - якобіан від F(X)
    global fun
    N = length(X); m = N-1;
    f = zeros(N,1);
    f(1) = (3 + 2.*X(1)).*X(1) - 2.*X(2) - 3;
    for i = 2:m
        f(i) = (3 + 2.*X(i)).*X(i) - X(i-1) - 2.*X(i+1) - 2;
    end
    f(N) = (3 + 2.*X(N)).*X(N) - X(m) - 4;
    if fun==0
        varargout = {f}; return
    end
    J = zeros(N);
    J(1,1) = 3 + 4.*X(1); J(1,2) = -2;
    for i = 2:m
        J(i,i-1) = -1; J(i,i)  = 3 + 4.*X(i); J(i,i+1) = -2;
    end
    J(N,m) = -1; J(N,N) = 3 + 4.*X(N);
    varargout = {f, J}; return
end

function print(txt, it, t, d, X, pxf)
%  друк
    global fun
    fprintf('%-18s %8u %8.4f %20.12e\n',txt,it,t,d);
    if pxf
        disp(X'); fun = 0; f = FdF(X); disp(f');
    end
end

function R = Dd(U,V)
% підфункція для обчислення матриці R,
%  де R := F(U;V) - узагальнені поділені різниці(I)
%     першого порядку від функції F(X) :
%            F(U;V)*(U - V) = F(U) - F(V).
    N = length(V); m = N-1;
    R = zeros(N);
    z3 =  2.*U(2) + 3; fl = (3 + 2.*U(1)).*U(1) - z3;
    z2 = (3 + 2.*V(1)).*V(1); fr = z2 - z3;
    R(1,1) = (fl - fr)./(U(1)-V(1));
    fl = fr; z3 = 2.*V(2) + 3; fr =  z2 - z3;
    R(1,2) = (fl - fr)./(U(2)-V(2));
    for i = 2:m
        im = i-1; ip = i+1;
        z1 = 2.*U(im); z2 = (3 + 2.*U(i)).*U(i); z3 = 2.*(U(ip) + 1);
        fl = -z1 + z2 - z3; z1 = 2.*V(im); fr = -z1 + z2 - z3;
        R(i,im) = (fl - fr)./(U(im)-V(im));
        fl = fr; z2 = (3 + 2.*V(i)).*V(i); fr = -z1 + z2 - z3;
        R(i,i)  = (fl - fr)./(U(i)-V(i));
        fl = fr; z3 = 2.*(V(ip) + 1); fr = -z1 + z2 - z3;
        R(i,ip) = (fl - fr)./(U(ip)-V(ip));
    end
    z2 = (3 + 2.*U(N)).*U(N); fl = z2 -4 - U(m);
    z3 = 4 + V(m); fr = z2 - z3;
    R(N,m) = (fl - fr)./(U(m)-V(m));
    fl = fr; z2 = (3 + 2.*V(N)).*V(N); fr = z2 - z3;
    R(N,N) = (fl - fr)./(U(N)-V(N));
    return
end