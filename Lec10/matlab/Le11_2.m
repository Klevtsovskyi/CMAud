%% Le11_2.m
% Розв'язання гіперболічного рівняння першого порядку
% з постійним коефіцієнтом   (неявна схема)
clear all, close all, clc

c = 1;
f  = @(x,t)(3+x.*(1-x)+t.*(3-2.*x));
u0 = @(x)(x.*(3-x));
g1 = @(t)(0);
a = 0; b = 3;

N = 301; x = linspace(a,b,N); h = x(2)-x(1); h2 = 0.5*h;
tau = 0.01; M = 500; tau2 = 0.5*tau;

C = c.*tau./h; B = 1./(1+C); A = C.*B; C = tau.*B;

t = 0; y = u0(x); U0 = y;
for j = 1:M
  t = t+tau; y(1) = g1(t);
  for k = 2:N
    y(k) = A.*y(k-1) + B.*y(k) + C.*f(x(k)-h2,t-tau2);
  end
  if j.*2 == M
     U1 = y;
  end
end

% побудова графіків розв'язку при t=0, M*tau/2, M*tau :
t=M.*tau; tau2=t./2; T1=num2str(tau2); T=num2str(t);
plot(x,U0,'b.-',x,U1,'k.--',x,y,'r*-');
xlabel('x'); ylabel('u');
title('Наближений розв''язок (неявна РС)');
legend('u(x,0)',['u(x,',T1,')'],['u(x,',T,')'],'Location','Best');
grid on; pause(); close all

u = @(x,t)((t+1).*x.*(3-x));
U = u(x,t); U0 = y;
load fileYJ y
plot(x, y, 'b.-', x, U, 'k.--');
xlabel('x'); ylabel('u');
title(['Набл.роз-к(явна РС) і точний розв''язок для t=',T]);
legend('явна','модель','Location','Best');
grid on; pause(); close all

plot(x, abs(y-U));
xlabel('x'); ylabel('u');
title(['Похибка |явна РС-модель| при t=',T]);
legend('|явна-модель|','Location','Best');
grid on; norm(y-U,inf); pause(); close all

plot(x, U0, 'b.-', x, U, 'k.--');
xlabel('x'); ylabel('u');
title(['Набл.роз-к(неявна РС) і точний розв''язок для t=',T]);
legend('неявна','модель','Location','Best');
grid on; pause(); close all

plot(x, abs(U0-U));
xlabel('x'); ylabel('u');
title(['Похибка |неявна РС-модель| при t=',T]);
legend('|неявна-модель|','Location','Best');
grid on; norm(U0-U,inf)