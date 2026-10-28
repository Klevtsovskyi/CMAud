%% Le11_3.m
% Розв'язання гіперболічного рівняння першого порядку
% зі зміним коефіцієнтом c(x,t)
clear all, close all, clc

c  = @(x,t)((t+0.5).*(0.6-x).^3);
f  = @(x,t)(sin(x-t));
u0 = @(x)(x.*(1-x));
g1 = @(t)(t);
g2 = @(t)(2.*t);
a = 0; b = 1;

N = 101; x = linspace(a,b,N); h = x(2)-x(1); N1 = N-1;
tau = 0.01; M = 200;

th = 0.5.*tau./h;

A = zeros(N,1); B = A; C = -ones(N,1); F = A;
t = 0; y = u0(x); U0 = y; U = zeros(N,M+1); U(:,1) = y;
for j = 1:M
  t = t+tau;
% формуєм коефіцієнти СЛАР   
  F(1) = -g1(t);
  for k = 2:N1
    z = c(x(k),t).*th; A(k) = -z; B(k) = z;
    F(k) = -(y(k)+tau.*f(x(k),t));
  end
  F(N) = -g2(t);
  
  [y, alfa] = m_progn(A,C,B,F);
  U(:,j+1) = y;
  
  if j.*2 == M
     U1 = y;
  end
end
clear alfa A C B F

% побудова графіків розв'язку при t=0, M*tau/2, M*tau :
t = M.*tau; z = t./2; T1 = num2str(z); T = num2str(t);
plot(x,U0,'b.-',x,U1,'k.--',x,y,'r*-');
xlabel('x'); ylabel('u');
title('Наближений розв''язок');
legend('u(x,0)',['u(x,',T1,')'],['u(x,',T,')'],'Location','Best');
grid on, pause(); close all, clear y U0 U1

tt = linspace(0,t,M+1);
[T,X] = meshgrid(tt,x);
surfc(X,T,U);
xlabel('x'); ylabel('t'); zlabel('u');
title('Наближений розв''язок u(X,T)');
shading flat
%colormap pink