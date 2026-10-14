function Pr08_1
clear all, close all, clc

F = @(x,U)([exp(-x)-U(1)]);

Y0 = [1];               % вектор початкових умов
X = linspace(0,5,501);  % рівномірна сітка з h=0.01
% розв'язання системи :
Y1 = Euler_ev(F, X, Y0);          % явний Ейлера
epsx = 1e-12; epsf = 1e-12; 
[Y2, efe] = Euler_nev(F, X, Y0, epsx, epsf);  % неявний Ейлера
if efe ~= 1
    efe
    return
end
[Y3, efs] = symm_nev(F, X, Y0, epsx, epsf);   % симетрична схема
if efs ~= 1
    efs
    return
end
u = @(x)(exp(-x).*(1+x));  Z = u(X);
plot(X, Y1, X, Y2, X, Y3, X, Z), grid on
title('u''(x)=exp(-x)-u'), xlabel('x'), ylabel('u'), 
legend('Ейлер яв','Ейлер неяв','симет неяв','точне','Location','Best')
pause(5), close all
Z1 = Z-Y1; Z2 = Z-Y2; Z3 = Z-Y3;
subplot(3,1,1); plot(X, Z1); grid on; title('u(x)-Y1'); hold on;
subplot(3,1,2); plot(X, Z2); grid on; title('u(x)-Y2'); hold on;
subplot(3,1,3); plot(X, Z3); grid on; title('u(x)-Y3'); hold off;
h = X(2)-X(1); hh = h.*h;
d1 = norm(Z1,inf); d2 = norm(Z2,inf); d3 = norm(Z3,inf);
fr = '%12s %12s %12s %12s %12s\n';
fprintf(fr,'h','h^2','||u-Y1||c','||u-Y2||c','||u-Y3||c');
fr = '%12.4e %12.4e %12.4e %12.4e %12.4e\n';
fprintf(fr,h,hh,d1,d2,d3);
end