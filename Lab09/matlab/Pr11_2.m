function Pr11_2
% п.2  K,p-неперервні, f-розривна :
% Модельна задача: при
%  K(x)=1, p(x)=0,
%  f(x)=-iif(x=0,0,100*sign(x)+exp(-x)),
%  a=-1, b=1; U(a)=0, U(b)=1;
% розв'язок буде
% U(x)=0.5*(1-99x)+50x^2*sign(x)+exp(x)
%      -0.5*((1+x)*exp(1)+(1-x)*exp(-1))
clear all, close all, clc
global X n h h05 h25 Trf nf
% Вхідні дані :
Kf = @(j)(1);
pf = @(j)(0);
a = -1; b = 1;
D = [0,1,0,0,1,1];
% Наближений розв'язок :
n = 201; X = linspace(a, b, n);
h = X(2)-X(1); h05 = 0.5*h; h25 = 0.5*h05;
% опис розміщення точки розриву f(x) на сітці X :
Trf = [101]; nf = 1;
Y = m_iim(Kf, pf, @ff, D);
% Точний розв'язок :
u = @(x)(0.5*(1-99*x)+50*x.*x.*sign(x)+exp(x)...
         -0.5*((1+x).*exp(1)+(1-x).*exp(-1)));
U = u(X);
t = 'K,p-неперервні,f-розривна';
view_solu(Y, U, 7, t, 5);
end

function e = ff(j)
% Функція f(x)=-iif(x=0,0,100*sign(x)+exp(x))
   global X n h25 Trf
   x = X(j); k = find(j == Trf(:));
   if isempty(k)
      %обчислення f(x) в неперервній точці по формулам (4.32)
      if j == 1
        x = x + h25;
      elseif j == n
        x = x - h25;
      end
      e = -100*sign(x) - exp(x);
   elseif k == 1
      %обчислення f(x) в точці 1-го розриву по формулі (4.33)
      e = -0.5*(exp(x-h25) + exp(x+h25));
   end
end