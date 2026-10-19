function Pr10_1
% Розв'язання модельної задачі для КЗ :
%   U"(x) - p(x)U(x) = -f(x),
%    p(x) >= 0,   x Є (a,b);
% D(1)*U'(a) + D(2)*U(a) = -D(3),
% D(4)*U'(b) + D(5)*U(b) = D(6).
% Для p(x)=4+x, f(x)=6+x(5+x)-(8+x)cos(x)^2, xє[0,3] i
%    D=[2,-1,-2,3,2,2,8.201583218752855]
%  розв'язок буде u(x)=x+sin(x)^2.
clear all, close all, clc
% Вхідні дані :
a = 0; b = 3;
p = @(x)(4 + x);
f = @(x)(6 + x.*(5 + x) - (8 + x).*cos(x).^2);
D = [0, 1, 0,...
     0, 1, 3+sin(b).^2];     % u(a)=-D(3), u(b)= D(6)
D = [1, 0,-1,...
     1, 0, 1+sin(2*6)];      % u'(a)=-D(3), u'(b)= D(6)
D = [2,-1, -2,...
     3, 2, 8.201583218752855];     
u = @(x)(x + sin(x).^2);   % точний розв'язок
fprintf('  h       ||y-u(X)||c\n');
fm = '%7.4f %14.6e\n'; k = 5;
for n = [31,301,3001,30001] 
    % Знаходження наближеного розв'язку :
    [X, Y] = m_grid2(p, f, a, b, D, n);
    h = X(2)-X(1);
    U = u(X);                  % проекція u(x) на сітку X
    m = min(U); M = max(U); I = [1:k:n]; k =10*k;
    close all;
    plot(X, Y, 'b -', X(I), U(I), 'k*'); grid on;
    s = sprintf('МСР O(h^2) h=%6.4f', h); title(s);
    de = 1e-1; axis([a-de, b+de, m-de, M+de]);
    xlabel('x'); ylabel('y,u(x)');
    legend('y','u(x)','Location','Best');
    pause(5); close all;
    % Візуалізація абсолютної похибки :
    U = abs(Y - U); M = max(U);
    plot(X, U, 'r', 'LineWidth', 1.5 ); grid on;
    title('Абсолютна похибка');
    xlabel('x'); ylabel('|y-u(x)|');
    fprintf(fm, h, M);
    pause(3);
end
end