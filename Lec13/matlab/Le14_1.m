function Le14_1
clear all, close all, clc
% Задача Діріхле для р-я Пуассона в квадратній області.
% Порівняння ітераційних процесів - 1.

a = 1; b = 1;
u = @(x,y)(1+15.*(x.*(a-x).*y.*(b-y)));
f = @(x,y)(30.*(x.*(a-x)+y.*(b-y)));

ep = 1e-8; itmax = 30000; Lep = log10(1./ep);
fm = '%12s %5u %5u%9.3f %13.5e';

N = [5; 11; 21; 41; 81];
nN = length(N); Kit = zeros(nN,3);
for n = 1 : nN
Nx = N(n); Ny = Nx; % кількість вузлів сітки
X = linspace(0, a, Nx); Y = linspace(0, b, Ny);
Hx = X(2)-X(1); Hy = Y(2)-Y(1); n0e = 2.*Lep./(pi.*Hx);

UU = zeros(Ny,Nx);
for j = 1 : Ny
   for i = 1 : Nx
      UU(j,i) = u(X(i),Y(j));
   end
end

mc = 'м.Якобі';
t1 = clock();
[U, it] = im_JZr(0, f, @muh, @muv, X, Y, 0, ep, itmax);
t = clock(); t = etime(t,t1);  Kit(n,1) = it;
fprintf('Nx=Ny=%4u Hx=Hy=%7.5f\n', Nx, Hx);
fprintf('Ітерац.метод Ітер. ІтерT Час(сек) max(max|U-u|) ro\n');
t1 = max(max(abs(U-UU))); n0 = fix(n0e./(pi.*Hx));
fprintf([fm,'\n'],mc, it, n0, t, t1);
viewU(2, X, Y, U, 15, UU);

mc = 'м.Зейделя';
t1 = clock();
[U, it] = im_JZr(1, f, @muh, @muv, X, Y, 1, ep, itmax);
t = clock(); t = etime(t,t1);  Kit(n,2) = it;
t1 = max(max(abs(U-UU))); n0 = fix(n0e./(2.*pi.*Hx));
fprintf([fm,'\n'],mc, it, n0, t, t1);
viewU(2, X, Y, U, 15, UU);

mc = 'м.релаксації';
ro = 1./(0.5+sin(0.5.*pi.*Hx));  % опт.параметр
t1 = clock();
[U, it] = im_JZr(2, f, @muh, @muv, X, Y, ro, ep, itmax);
t = clock(); t = etime(t,t1);  Kit(n,3) = it;
t1 = max(max(abs(U-UU))); n0 = fix(n0e);
fprintf([fm,' %7.5f\n'], mc, it, n0, t, t1, ro);
viewU(2, X, Y, U, 15, UU);
end

pause, close all
clear U UU
N = N - 1; M = linspace(N(1),N(end),101);
r = 2.*Lep./pi.*M; J = r./pi.*M; Z = 0.5.*J;
plot(N,Kit(:,1),'ks ','LineWidth',2); hold on;
plot(N,Kit(:,2),'ko ','LineWidth',2); hold on;
plot(N,Kit(:,3),'b* ','LineWidth',2); hold on;
plot(M,J,'r -.','LineWidth',2); hold on;
plot(M,Z,'r -','LineWidth',2); hold on;
plot(M,r,'r --','LineWidth',2); hold off;
legend('Якобі','Зейделя','релаксації',...
       'Як n0','Зе n0','ре n0','Location','Best');
grid on; xlabel('h^{-1}'); ylabel('Kit');
title('Число ітерацій від h');

end

function [g1, g2] = muh(x)
   g1 = 1.0; g2 = 1.0;
end

function [g1, g2] = muv(y)
   g1 = 1.0; g2 = 1.0;
end
