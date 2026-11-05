function Pr12_1
clear all, close all, clc
% п.1 модельна задача
% параметри задачі :
a = 0; b = 1;
p = @(x,t)((2+x).*t);
f = @(x,t)((2+x.*(x-2)+(4-x.*x).*x.*t).*exp(-t));
u0 = @(x)(x.*(2-x));
v = [1, 1];
% параметри різницевої схеми :
n = 101; X = linspace(a, b, n);
tau = 4e-5; M = 78000; J = [0 : 1200 : M];
% розв'язання задачі :
tim1 = clock();
U =  fdm_P10(p, f, u0, @Ga, @Gb, v, X, J, tau);
tim2 = clock(); C = etime(tim2,tim1);
fprintf('Тривалість розрахунку =%13.3f сек\n',C);
J = tau.*J; m = length(J); T = J(m);
% точний розв'язок :
u = @(x,t)(x.*(2-x).*exp(-t));
Uu = U;
for i = 1 : n
   x = X(i);
   for k = 0 : m-1
      j = k+1; t = J(j); Uu(j,i) = u(x,t);
   end
end
[Mx, Mt] = meshgrid(X, J); % побудова сіток для 3D-графіків
% виведення графіків: розв'язку U і точного u(x,t)
de = [a, b, 0, T, -0.2, 1.2];
mesh(Mx, Mt, U); axis(de);
shading interp, colormap cool;
xlabel('x'), ylabel('t'), zlabel('U')
title('U - наближений');
pause(), close all
mesh(Mx, Mt, Uu); axis(de);
shading interp, colormap cool;
xlabel('x'); ylabel('t'); zlabel('u');
title('u(x,t) - точний');
pause(), close all
% виведення графіків похибки U-Uu
Uu = U-Uu;
C = norm(Uu,inf); 
fprintf('Похибка ||U-u||c=%13.5e\n',C);
C = max(max(abs(Uu))); 
fprintf('max(max(|U-u|)) =%13.5e\n',C);
mesh(Mx, Mt, Uu);
shading interp, colormap jet;
xlabel('x'); ylabel('t'); zlabel('U-u');
title('Похибка U-u(x,t)');
end

function gr = Ga(t)
% Опис [g(t), мю(t)] для x=a
% КУ I   роду  va=0 :                  1*u(a,t) = мю(t)
% КУ II  роду  va=1 : -du(a,t)/dx +    0*u(a,t) = мю(t)
% КУ III роду  va=1 : -du(a,t)/dx + g(t)*u(a,t) = мю(t)
  gr = [ 1, -2*exp(-t)];          % III р
 %gr = [ 0, -2*exp(-t)];          % II  р
 %gr = [ 1, 0];                   % I   р
end

function qr = Gb(t)
% Опис [q(t), ню(t)] для x=b
% КУ I   роду  vb=0 :                  1*u(b,t) = ню(t)
% КУ II  роду  vb=1 :  du(b,t)/dx +    0*u(b,t) = ню(t)
% КУ III роду  vb=1 :  du(b,t)/dx + q(t)*u(b,t) = ню(t)
   qr = [ 2, 2.*exp(-t)];         % III р
  %qr = [ 1, exp(-t)];            % II  р
  %qr = [ 1, 0 ];                 % I   р
end