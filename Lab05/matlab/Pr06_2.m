% Pr06_2.m
clear all, close all, clc, warning off

a = 2./23; a1 = 1 - a.*a;
f = @(x)(sqrt(1-a1.*sin(x).^2));
tol = 1e-10;
A = 0; B = 0.5.*pi;

N = [2,4,6,7];
E = 4.*quadl(f, A, B, tol);
fprintf('E= %.11e%s\n',E,' - вбудована м.Лобатто');
fprintf(' n КФ Гаусса	    КФ Чебишева	     КФ Маркова\n');
for i = 1 : 4
   n = N(i);
   g = 4.*i_gauss(A, B, f, n);
   c = 4.*i_cheb(A, B, f, n);
   m = 4.*i_markov(A, B, f, n);
   fprintf('%2d %.9e %.9e %.9e\n', n, g, c, m);
end