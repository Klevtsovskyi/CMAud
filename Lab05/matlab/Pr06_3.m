% Pr06_3.m
clear all, close all, clc, warning off

a = -1; b = 2;
f = @(x)(x.*abs(x));
fm= @(x)(-x.*x);
fp= @(x)(x.*x);

S = 7/3; Nn = 512; Nm =5; Np = 11;
fprintf('            S точне= %.9e\n', S);
Ss = i_simpn(a, 0, fm, Nm)+i_simpn(0, b, fp, Np);
fprintf('          S вар. A)= %.9e\n', Ss);
fprintf('%s\n', 'вар. B)');
fprintf('  n КФ трапецій      КФ центр.прямок.\n');
n = 1;
while n < Nn
   n = 2*n;
   t = i_trapn(a, b, f, n);
   z = i_crectn(a, b, f, n);
   fprintf('%3d %.9e %.9e\n', n, t, z);
end
