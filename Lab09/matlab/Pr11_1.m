function Pr11_1
% п.1  Вхідні дані неперервні :
% Модельна задача: при
%  K(x)=1+x, p(x)=exp(-x),
%  f(x)=sin(x)*(1+x+exp(-x))-cos(x),
%  a=1, b=3;
%  U(a)=sin(a), U(b)=sin(b);
% розв'язок буде U(x)=sin(x).  
clear all, close all, clc
global X n h h05 h25
% Вхідні дані :
a = 1; b = 3;
D=[0,1,sin(a),0,1,sin(b)];                                  %  Iр,  Iр 
D=[1,0,-(1+a).*cos(a),0,1,sin(b)];                          % IIр,  Iр
D=[1,2,-(1+a).*cos(a)+2*sin(a),0,1,sin(b)];                 %IIIр,  Iр
D=[0,1,sin(a),1,0,(1+b).*cos(b)];                           %  Iр, IIр
D=[1,0,-(1+a).*cos(a),1,0,(1+b).*cos(b)];                   % IIр, IIр
D=[1,2,-(1+a).*cos(a)+2*sin(a),1,0,(1+b).*cos(b)];          %IIIр, IIр
D=[0,1,sin(a),1,3,(1+b).*cos(b)+3*sin(b)];                  %  Iр,IIIр
D=[1,0,-(1+a).*cos(a),1,3,(1+b).*cos(b)+3*sin(b)];          % IIр,IIIр
D=[1,2,-(1+a).*cos(a)+2*sin(a),1,3,(1+b).*cos(b)+3*sin(b)]; %IIIр,IIIр
% Наближений розв'язок :
n = 201; X = linspace(a, b, n);
h = X(2)-X(1); h05 = 0.5*h; h25 = 0.5*h05;
Y = m_iim(@Kf, @pf, @ff, D);
% Точний розв'язок :
u = @(x)(sin(x)); U = u(X);
t = 'K,p,f-неперервні';
view_solu(Y, U, 7, t, 5);
end

function e = Kf(j)
% Функція K(x)=1+x
    global X h05
    x = X(j) - h05;
    e = 1 + x;
end

function e = pf(j)
% Функція p(x)=exp(-x)
    global X n h25
    x = X(j);
    if j == 1
        x = x + h25;
    elseif j == n
        x = x - h25;
    end
    e = exp(-x);
end

function e = ff(j)
% Функція f(x)=sin(x)*(1+x+exp(-x))-cos(x)
    global X n h25
    x = X(j);
    if j == 1
        x = x + h25;
    elseif j == n
        x = x - h25;
    end
    e = sin(x).*(1 + x + exp(-x)) - cos(x);
end