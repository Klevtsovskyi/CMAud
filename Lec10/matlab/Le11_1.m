%% Le11_1.m
% Розв'язання гіперболічного рівняння першого порядку
% з постійним коефіцієнтом   (явна схема)
clear all, close all, clc

c = 1;
f  = @(x,t)(3+x.*(1-x)+t.*(3-2.*x));
u0 = @(x)(x.*(3-x));
g1 = @(t)(0);
a = 0; b = 3;

N = 301; x = linspace(a,b,N); h = x(2)-x(1); h2 = 0.5*h;
tau = 0.001; M = 5000; tau2 = 0.5*tau;

A = c.*tau./h; B = 1-A;

t = 0; y = u0(x); U0 = y;
for j = 1:M
  t = t+tau;
  yr = y(1); y(1) = g1(t);
  for k = 2:N
    yl = yr; yr = y(k);
    y(k) = A.*yl + B.*yr + tau.*f(x(k)-h2,t-tau2);
  end
  if j.*2 == M
     U1 = y;
  end
end
save fileYJ y

% побудова графіків розв'язку при t=0, M*tau/2, M*tau :
t=M.*tau; tau2=t./2; T1=num2str(tau2); T=num2str(t);
plot(x,U0,'b.-',x,U1,'k.--',x,y,'r*-');
xlabel('x'); ylabel('u');
title('Наближений розв''язок (явна РС)');
legend('u(x,0)',['u(x,',T1,')'],['u(x,',T,')'],'Location','Best');
grid on