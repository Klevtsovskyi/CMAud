% Pr05_1.m
clear all, close all, clc, warning off

u   = @(x)(sin(x));
du  = @(x)(cos(x));
d2u = @(x)(-sin(x));

X = linspace(0,3,31); n = length(X); n1 = n-1;
h = X(2)-X(1); hh = 2*h; h2 = h.*h;
U = u(X); U1 = du(X); U2 = d2u(X);
Yl=0.*X; Yr=Yl; Yz=Yl; Ylr=Yl; R=Yl;
for i=2:n1
    im = i-1; ip = i+1;
    Yl(i) = (U(i)-U(im))./h;           % ліва
    Yr(i) = (U(ip)-U(i))./h;           % права
    Yz(i) = (U(ip)-U(im))./hh;         % центральна
    Ylr(i) = (U(im)-2*U(i)+U(ip))./h2; % друга різницева
end
Yl(n) = (U(n)-U(n1))./h;
Yr(1) = (U(2)-U(1))./h;

subplot(1,2,1);
I=2:n;
plot(X(I),U1(I),'k :',X(I),Yl(I),'b');
grid on; xlabel('X'); ylabel('U''');
title('U''(X) і ліва різницева');
subplot(1,2,2);
frm1 = '%-10s %12s %12s\n'; frm2 = '%10.8f %12.10f %12.10f\n';
fprintf(frm1,'Крок h','Похибка теор','Похибка розр');
R(I) = U1(I)-Yl(I); r = norm(R(I),inf); fprintf(frm2,h,h/2,r);
plot(X(I),R(I));
grid on; xlabel('X'); ylabel('R(X)');
title('Похибка R(X)'); pause; close all

subplot(1,2,1);
I=1:n1;
plot(X(I),U1(I),'k :',X(I),Yr(I),'b');
grid on; xlabel('X'); ylabel('U''');
title('U''(X) і права різницева');
subplot(1,2,2);
R(I) = U1(I)-Yr(I); r = norm(R(I),inf); fprintf(frm2,h,h/2,r);
plot(X(I),R(I));
grid on; xlabel('X'); ylabel('R(X)');
title('Похибка R(X)'); pause; close all

subplot(1,2,1);
I=2:n1;
plot(X(I),U1(I),'k :',X(I),Yz(I),'b');
grid on; xlabel('X'); ylabel('U''');
title('U''(X) і центральна');
subplot(1,2,2);
R(I) = U1(I)-Yz(I); r = norm(R(I),inf); fprintf(frm2,h,h2/6,r);
plot(X(I),R(I));
grid on; xlabel('X'); ylabel('R(X)');
title('Похибка R(X)'); pause; close all

subplot(1,2,1);
plot(X(I),U2(I),'k :',X(I),Ylr(I),'b');
grid on; xlabel('X'); ylabel('U"');
title('U"(X) і друга різницева');
subplot(1,2,2);
R(I) = U2(I)-Ylr(I); r = norm(R(I),inf); fprintf(frm2,h,h2/12,r);
plot(X(I),R(I));
grid on; xlabel('X'); ylabel('R(X)');
title('Похибка R(X)'); pause; close all

P = polyfit(X,U,n1);
dP = polyder(P);
dY = polyval(dP,X);
d2P = polyder(dP);
d2Y = polyval(d2P,X);

subplot(1,2,1);
plot(X,U1,'k :',X,dY,'b');
grid on; xlabel('X'); ylabel('U''');
title('U''(X) і dP');
subplot(1,2,2);
R = U1-dY; r = norm(R,inf);
plot(X,R);
grid on; xlabel('X'); ylabel('R(X)');
title('Похибка R(X)'); pause; close all

subplot(1,2,1);
plot(X,U2,'k :',X,d2Y,'b');
grid on; xlabel('X'); ylabel('U"');
title('U"(X) і d2P');
subplot(1,2,2);
R = U2-d2Y; r = norm(R,inf);
plot(X,R);
grid on; xlabel('X'); ylabel('R(X)');
title('Похибка R(X)'); pause; close all