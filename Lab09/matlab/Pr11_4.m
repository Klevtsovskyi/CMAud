function Pr11_4
% п.4  p,f-неперервні, K-розрив 1-го роду :
% Крайова залача :
%  (K(x)*u'(x))' = 0, xє(0,1)
%  K(x)=K1, 0<x<c або =K2, c<x<1; c = X(m)+0.5h, 1<m<n;
%  U(0)=1, U(1)=0.
clear all, close all, clc
global X n h h05 h25 TrK nK K1 K2 c al be
% Вхідні дані :
pf = @(j)(0);
ff = @(j)(0);
K1 = 2; K2 = 3.5; a = 0; b = 1; D = [0,1,1,0,1,0];
% Наближений розв'язок :
n = 101; X = linspace(a, b, n);
h = X(2)-X(1); h05 = 0.5*h; h25 = 0.5*h05;
% опис точки розриву K(x) на сітці X :
TrK = [39]; nK = 1; c = X(TrK(1)) + h05; t = K1./K2;
al = 1/(t + (1-t).*c); be = al.*t;
Y = m_iim(@Kf, pf, ff, D);
U = X;
for k = 1:n
   U(k) = u(X(k));  % точний розв'язок
end
t = 'p,f-неперервні,K-розривна';
view_solu(Y, U, 21, t, 5);
end

function e = Kf(j)
% Функція K(x)=K1, 0<x<c,
%             =K2, c<x<1; c = X(m)+h05, 1<m<n
   global TrK K1 K2
   m = TrK(1)+1;
   if j < m
      e = K1;
   elseif j == m
      e = 0.5*(K1 + K2); 
   else
      e = K2;
   end
end

function e = u(x)
% Функція u(x)
    global c al be
    if x <= c
        e = 1 - al.*x;
    else
        e = be.*(1 - x);
    end
end