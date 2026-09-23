% Le04_1.m
% Приклади числового диференціювання
clear all, close all, clc, warning off

n = 101; x = linspace(0, 2*pi, n);
n1 = n - 1; h = x(2) - x(1); hh = h.*h; h2 = 2*h;
y = sin(x); dy = cos(x);
% 1-ша похідна
Dy = diff(y);           
%dydx = diff(y)./diff(x);        %  у випадку нерівномірної сітки  
dydx = Dy./h;                    %  для точок x(1:n1) - права РП
I = 1:n1;
dydxc = (y(3:n) - y(1:n-2))./h2; %  для точок x(2:n1) - центральна РП
subplot(3,1,1);
plot(x(I),dydx - dy(I)); grid on; hold on;
title('Похибка правої РП');
subplot(3,1,2); I = 2:n1;
plot(x(I),dydxc - dy(I)); grid on; hold on;
title('Похибка центральної РП');
% 2-га похідна
%  для точок x(2:n1) - 2-га РП
d2ydx2 = zeros(1,n-2);
for k = 2 : n1
   k1 = k - 1;
   d2ydx2(k1) = (y(k1)-2*y(k)+y(k+1))./hh;
end
subplot(3,1,3);
plot(x(I),d2ydx2 + y(I)); grid on; hold off;
title('Похибка 2-ї РП')
pause; close all; clear x y dydx dydxc d2ydx2

f = @(x,y)(16 - x.^4 - y.^4);
x = linspace(-2, 2, 41); hx = x(2) - x(1);
y = linspace(-2, 2, 21); hy = y(2) - y(1);
[X, Y] = meshgrid(x, y); Z = f(X, Y);
[Px, Py] = gradient(Z, hx, hy);
contour(Z, 10); hold on;
quiver(Px, Py); hold off;
s = char(f); title(s);
pause; close all; clear Px Py

surfc(X, Y, Z); grid on; title(s)
pause; close all;
L = 4.*del2(Z, hx, hy);
surfc(X, Y, L); grid on;
s = ['{\Delta}',s]; title(s);