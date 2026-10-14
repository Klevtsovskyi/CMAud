function Pr08_2
clear all, close all, clc

F = @(x,U)([exp(-x)-U(1)]);

Y0 = [1];               % вектор початкових умов
X = linspace(0,5,501);  % рівномірна сітка з h=0.01
Y1 = rk2_ev(F, X, Y0, 1);           % 2-стадійний Р-К на рівн.сітці
[T1, Y2] = ode23(F, X, Y0); Y2=Y2'; % солвер
u = @(x)(exp(-x).*(1+x)); Z = u(X); % точний
plot(X, Y1, X, Y2, X, Z), grid on
title('u''(x)=exp(-x)-u'), xlabel('x'), ylabel('u'), 
legend('rk2','ode23','точне','Location','Best')
pause(5), close all
Z1 = Z-Y1; Z2 = Z-Y2;
subplot(2,1,1); plot(X, Z1); grid on; title('u(x)-Y1'); hold on;
subplot(2,1,2); plot(X, Z2); grid on; title('u(x)-Y2'); hold off;
h = X(2)-X(1); h2 = h.*h; h4 = h2.*h2;
d1 = norm(Z1,inf); d2 = norm(Z2,inf);
fr1 = '%12s %12s %12s\n';
fprintf(fr1,'h^2','||u-Y1||c','||u-Y2||c');
fr2 = '%12.4e %12.4e %12.4e\n';
fprintf(fr2,h2,d1,d2);
pause(5), close all

Y1 = rk4_ev(F, X, Y0);              % 4-стадійний Р-К на рівн.сітці
[T1, Y2] = ode45(F, X, Y0); Y2=Y2'; % солвер
plot(X, Y1, X, Y2, X, Z), grid on
title('u''(x)=exp(-x)-u'), xlabel('x'), ylabel('u'), 
legend('rk4','ode45','точне','Location','Best')
pause(5), close all
Z1 = Z-Y1; Z2 = Z-Y2;
subplot(2,1,1); plot(X, Z1); grid on; title('u(x)-Y1'); hold on;
subplot(2,1,2); plot(X, Z2); grid on; title('u(x)-Y2'); hold off;
d1 = norm(Z1,inf); d2 = norm(Z2,inf);
fprintf(fr1,'h^4','||u-Y1||c','||u-Y2||c');
fprintf(fr2,h4,d1,d2);
end