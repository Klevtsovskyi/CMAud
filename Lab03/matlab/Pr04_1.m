% Pr04_1.m
clear all, close all, clc, warning off

u   = @(t)(exp(sin(pi.*t))./t.^3);
      
n = 11; nm = n-1; s = 5; N = nm*s+1;
a = 1.5; b = 3.5;
Xr = linspace(a, b, n);        % сітка рівномірна  
                               % сітка нерівномірна :
Xn = sort(Xr(1)+(Xr(n)-Xr(1)).*rand(1,n));
Ur = u(Xr); Un = u(Xn);        % значення u(t) на сітках
x = linspace(a, b, N);         % сітка рівномірна для обчислень "згущена"
I = 1:s:N;                     % індексний масив
Pr = polyfit(Xr, Ur, nm);      % побудова інтерп.полінома Pr на Xr
Yr = polyval(Pr, x);           % значення ІП Pr в точках x
Pn = polyfit(Xn, Un, nm);      % побудова інтерп.полінома Pn на Xn
Yn = polyval(Pn, x);           % значення ІП Pn в точках x

subplot(1,2,1);
plot(Xr,Ur,'ko :',x,Yr,'b');
grid on; xlabel('Xr,x'); ylabel('Ur,Yr');
title('сітка рівномірна');
subplot(1,2,2);
Rr = Yr(I) - Ur;               % похибка 
plot(Xr, Rr);
grid on; xlabel('Xr'); ylabel('Rr');
title('Похибка Rr=Yr-Ur'); pause(5); close all

subplot(1,2,1);
plot(Xn,Un,'ko :',x,Yn,'b');
grid on; xlabel('Xn,x'); ylabel('Un,Yn');
title('сітка нерівномірна');
subplot(1,2,2);
Rn = Yn(I) - Un;               % похибка 
plot(Xn, Rn);
grid on; xlabel('Xn'); ylabel('Rn');
title('Похибка Rn=Yn-Un'); pause(5); close all

plot(Xr,Rr,'k',Xn,Rn,'b');
grid on; xlabel('Xr,Xn'); ylabel('Rr,Rn');
title('похибки для різних сіток');