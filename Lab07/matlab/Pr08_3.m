% Pr08_3.m
clear all, close all, clc, warning off

% Розв'язання задачі Коші
%   u"+3u'+7u=sin(u);
%    u(0)=3,u'(0)=-1.
f1 = @(u1,u2)(u2);                  % опис функції f1(u1,u2)
f2 = @(u1,u2)(sin(u1)-3.*u2-7.*u1); % f2(u1,u2)
f = @(t,U)[f1(U(1),U(2)); f2(U(1),U(2))];       % f(t,U)
U0 = [3;-1];                        % вектор початкових умов
% Розв'язання задачі Коші для системи ЗДР
[T, U] = ode45(f, [0 10], U0);
u = U(:,1); du = U(:,2);
clear U U0

% Побудова графіків u(t), u'(t)
subplot(2,2,1);
plot(T, u, 'b', 'LineWidth', 2), hold on
plot(T, du, 'r'), hold off, grid on
s = 'u(t)'; sn = 'u''(t)';
legend(s, sn, 'Location', 'Best')
title('4u"+5u''+9u=cos(u); u(0)=4, u''(0)=-3')
xlabel('t'); ylabel([s,',',sn]); clear T;

s = 'Фазова траекторія';
subplot(2,2,2);
plot(u, du, 'k', 'LineWidth', 2);
axis equal; grid on;
title(s); xlabel('u'), ylabel('u''');

sn = 'Поле напрямків';
subplot(2,2,3);
m1 = min(u); M1 = max(u);
x = linspace(m1,M1,11);
m2 = min(du); M2 = max(du);
y = linspace(m2,M2,11);
[X, Y] = meshgrid(x, y);
F1 = f1(X, Y); F2 = f2(X, Y);
clear x y m1 m2 M1 M2;

quiver(X, Y, F1, F2)
axis equal; title(sn);

s = [s,'+',sn];
subplot(2,2,4);
quiver(X, Y, F1, F2)
axis equal; hold on
plot(u, du, 'k', 'LineWidth', 2);
axis equal; hold off; grid on;
title(s); xlabel('u'), ylabel('u''')