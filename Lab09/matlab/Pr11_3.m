function Pr11_3
% п.3  K,f-неперервні, p-розривна :
% Крайова залача :
%  u"(x)-p(x)*u(x)=-sin(pi*x), xє(0,1)
%  p(x)=p1, 0<=x<c,
%      =p2, c<x<=1; c = X(m), 1<m<n;
%  U(0)=1, U'(1)+3*U(1)=4.
clear all, close all, clc
global X n h h05 h25 Trp np p1 p2
% Вхідні дані :
Kf = @(j)(1);
p1 = 2; p2 = 10;
a = 0; b = 1;
D = [0,1,1,1,3,4];
% Наближений розв'язок :
n = 101; X = linspace(a, b, n);
h = X(2)-X(1); h05 = 0.5*h; h25 = 0.5*h05;
% опис розміщення точки розриву p(x) на сітці X :
Trp = [31]; np = 1;
Y = m_iim(Kf, @pf, @ff, D);
t = 'K,f-неперервні,p-розривна';
view_sol(Y, t);
end

function e = pf(j)
% Функція p(x)=p1, 0 <=x < c,
%             =p2, c < x <=1; c = X(m), 1<m<n
   global X n h25 Trp p1 p2
   x = X(j); k = find(j == Trp(:));
   if isempty(k)
      %обчислення p(x) в неперервній точці по формулам (4.32)
      if j == 1
        x = x + h25;
      elseif j == n
        x = x - h25;
      end
      if j < Trp(1)
         e = p1;
      else
         e = p2;
      end
   elseif k == 1
      %обчислення p(x) в точці 1-го розриву по формулі (4.33)
      e = 0.5*(p1 + p2);
   end
end

function e = ff(j)
% Функція f(x)=sin(pi*x)
    global X n h25
    x = X(j);
    if j == 1
        x = x + h25;
    elseif j == n
        x = x - h25;
    end
    e = sin(pi.*x);
end