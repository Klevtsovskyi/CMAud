function Pr13_3
clear all, close all, clc, warning off
% п.3 модельна задача
% параметри задачі :
a = 0; b = 1;
K = @(u)(1.0);
f = @(u)((1+u.*u).*(1-2.*u));
u0 = @(x)(tan(x));
% параметри різницевої схеми :
n = 201; X = linspace(a, b, n);
tau = 1e-3; M = 400; J = [0 : 10 : M];
% розв'язання задачі :
tim1 = clock();
Y =  fdm_P1n3(K, f, u0, @mn, X, J, tau);
tim2 = clock(); C = etime(tim2,tim1);
fprintf('Тривалість розрахунку =%13.3f сек\n',C);
J = tau.*J; m = length(J); T = J(m);
% точний розв'язок :
u = @(x,t)(tan(x+t)); U = Y;
for i = 1 : n
   x = X(i);
   for j = 1 : m
      U(j,i) = u(x,J(j));
   end
end
ti1 = 'Y - наближений'; z1 = 'Y';
ti2 = 'U - точний'; z2 = 'U';
tp = 'Похибка'; z3 = 'Y-U'; ti3 = [tp,' ',z3];
view_solu(X, J, Y, U, ti1, z1, ti2, z2, tp, ti3, z3, 5);
end

function e = mn(t)
% Опис [мю(t), ню(t)] для x=a, x=b
   e = [ tan(t), tan(1+t)];
end
