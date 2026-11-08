function Le13_1
clear all, close all, clc
% ППС для модельної задачі
% параметри задачі :
a1 = 0; b1 = 1; a2 = 0; b2 = 2;
f  = @(x1,x2,t)(exp(-t).*(4+x1.*x1+x2.*x2)-7); % права частина ДР
u0 = @(x1,x2)(x1.*x1+x2.*x2);                  % початкова умова КЗ
% параметри різницевої схеми :
N1 = 21; X1 = linspace(a1, b1, N1);
N2 = 41; X2 = linspace(a2, b2, N2);
tau = 0.0025; M = 2000; J = [];
% розв'язання задачі :
tim1 = clock();
Y = fmd_P2PR(f, @muh, @muv, u0, X1, X2, tau, M, J);
tim2 = clock(); C = etime(tim2,tim1);
fprintf('Тривалість розрахунку =%13.3f сек\n',C);
% точний розв'язок :
u = @(x1,x2,t)(t+(x1.*x1+x2.*x2).*(2-exp(-t)));
U = zeros(N2,N1); T = tau.*M;
for k1 = 1 : N1
   x1 = X1(k1);
   for k2 = 1 : N2
      U(k2,k1) = u(x1,X2(k2),T);
   end
end
ti1 = 'Y - наближений'; z1 = 'Y';
ti2 = 'U - точний'; z2 = 'U';
tp = 'Похибка'; z3 = 'Y-U'; ti3 = [tp,' ',z3];
view_solu(X1, X2, Y, z1, ti1, 0, 5, U, z2, ti2, tp, z3, ti3);
end

function [ga, gb] = muh(x1, t)
% Опис [мю(x1,a2,t), мю(x1,b2,t)]
   xx = x1.*x1; e = 2 - exp(-t);
   ga = t + xx.*e; gb = t + (xx + 4).*e;
end

function [ga, gb] = muv(x2, t)
% Опис [мю(a1,x2,t), мю(b1,x2,t)]
   xx = x2.*x2; e = 2 - exp(-t);
   ga = t + xx.*e; gb = t + (1 + xx).*e;
end
