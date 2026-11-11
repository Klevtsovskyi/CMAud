function Pr13_2
clear all, close all, clc, warning off
% п.2 модельна задача
% параметри задачі :
a = 0; b = 1;
fi = @(u)(u);
u0 = @(x)(tan(x));
% параметри різницевої схеми :
n = 201; X = linspace(a, b, n);
tau = 2.5e-5; M = 16000; J = [0 : 320 : M];
ep = 1e-10; sm = 10;
% розв'язання задачі :
tim1 = clock();
Y =  fdm_P1n2(fi, @fdf, u0, @mn, X, J, tau, ep, sm);
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

function e = fdf(u)
% Опис [f(u), d(f(u))/du, d(фі(u))/du]
    v = 1 + u.*u; z = 2.*u; w = 1 - z;
    e = [ v.*w, z.*w-2.*v, 1];
end

function e = mn(t)
% Опис [мю(t), ню(t)] для x=a, x=b
   e = [ tan(t), tan(1+t)];
end
