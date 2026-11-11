function Pr13_1
clear all, close all, clc, warning off
% п.1 модельна задача
% параметри задачі :
a = 0; b = 1;
K = @(x,t,u)(1.0);
u0 = @(x)(tan(x));
% параметри різницевої схеми :
n = 201; X = linspace(a, b, n);
tau = 2.5e-5; M = 16000; J = [0 : 320 : M];
ep = 1e-10; sm = 10;
% розв'язання задачі :
tim1 = clock();
Y =  fdm_P1n1(K, @cf, u0, @mn, X, J, tau, ep, sm);
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

function e = cf(x,t,u)
% Опис [c(x,t,u), f(x,t,u)]
    e = [ 1, (1+u.*u).*(1-2.*u)];
end

function e = mn(t)
% Опис [мю(t), ню(t)] для x=a, x=b
   e = [ tan(t), tan(1+t)];
end
