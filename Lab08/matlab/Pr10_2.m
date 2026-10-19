function Pr10_2
% Розв'язання модельної задачі для КЗ :
%   U"(x) - p*U(x) = -f(x),
% p є Const, p >= 0, x є (a,b);
% D(1)*U'(a) + D(2)*U(a) = -D(3),
% D(4)*U'(b) + D(5)*U(b) =  D(6).
% Для p=4, f(x)=6+4x-8cos(x)^2, xє[0,3] i
%    D=[2,-1,-2,3,2,8.201583218752855]
%  розв'язок буде u(x)=x+sin(x)^2.
clear all, close all, clc
% Вхідні дані :
a = 0; b = 3;
p = @(x)(4); P = 4;
f = @(x)(6+4*x-8*cos(x).^2);
u = @(x)(x + sin(x).^2);   % точний розв'язок
%D = [0, 1, 0, 0, 1, u(b)];              %  I р,  I р
%D = [0, 1, 0, 1, 0, 1+sin(2*b)];        %  I р, II р
%D = [0, 1, 0, 3, 2, 8.201583218752855]; %  I р,III р
%D = [1, 0,-1, 0, 1, u(b)];              % II р,  I р
%D = [1, 0,-1, 1, 0, 1+sin(2*b)];        % II р, II р
%D = [1, 0,-1, 3, 2, 8.201583218752855]; % II р,III р
%D = [2, 1,-2, 0, 1, u(b)];              %III р,  I р
%D = [2,-1,-2, 1, 0, 1+sin(2*b)];        %III р, II р
 D = [2,-1,-2, 3, 2, 8.201583218752855]; %III р,III р
n = 301; k = 50;
% Знаходження наближеного розв'язку :
[X, Y1] = m_grid2(p, f, a, b, D, n);
h = X(2)-X(1); U = u(X);   % проекція u(x) на сітку X
m = min(U); M = max(U); I = [1:k:n];
de = 1e-1; de = [a-de, b+de, m-de, M+de];
plot(X, Y1, 'b', X(I), U(I), 'k*', 'LineWidth', 2);
grid on; title('МСР O(h^2)'); axis(de);
xlabel('x'); ylabel('Y1,u(x)');
legend('Y1','u(x)','Location','Best');
pause(5); close all;
E1 = abs(Y1 - U); mE = max(E1);      % абсолютна похибка
fprintf('  h       ||Y-u(X)||c\n');
fm = '%7.4f %14.6e\n';
fprintf(fm, h, mE);

[X, Y2] = m_grid4(P, @fd2f, @df1, a, b, D, n);
plot(X, Y2, 'r', X(I), U(I), 'k*', 'LineWidth', 2);
grid on; title('МСР O(h^4)'); axis(de);
xlabel('x'); ylabel('Y2,u(x)');
legend('Y2','u(x)','Location','Best');
pause(5); close all;
% Візуалізація абсолютних похибок :
E2 = abs(Y2 - U); mE = max(E2);
fprintf(fm, h, mE);
plot(X, E1, 'b', 'LineWidth', 1.5 ); hold on;
plot(X, 1e5.*E2, 'r', 'LineWidth', 1.5 ); grid on;
fm = '|{\psi}|';
title(['Абсолютна похибка ',fm,':=|Y-u(x)|']); hold off;
xlabel('x'); ylabel(fm);
fm = [fm,' - O(h^'];
legend([fm,'2)'],['10^5*',fm,'4)'],'Location','Best');
format long e
disp('     U                         Y1                        Y2');
disp([U(I)',Y1(I)',Y2(I)']);
end

function [f, f2] = fd2f(x)
% Опис f(x), f"(x)
   f  = 6+4*x-8*cos(x).^2;
   f2 = 16*cos(2*x);
end

function df = df1(x)
% Опис f'(x)
   df = 4+8*sin(2*x);
end