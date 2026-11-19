function Pr14_1
clear all, close all, clc
%  модельна задача :
%   Розв'язування крайової задачі для рівняння Пуассона
%   в області складеної з класичних геометричних елементів
%   (прямокутників, трапецій, трикутників)

% параметри задачі і "опис" геометрії області :
a = 1; b = 1;
Nx = 201; Ny = 101; % кількість вузлів узгодженої сітки
Nx1 = Nx-1; Ny1 = Ny-1;
X = linspace(0, a, Nx); Y = linspace(0, b, Ny);
Hx = X(2)-X(1); Hy = Y(2)-Y(1);
Nx11 = 61; Nx21 = 101; Ny11 = 31;

% ітераційні параметри м.релаксації :
ro = 1.965; ep = 1e-8; itmax = 30000;
% допоміжні величини :
HxH = 1./(Hx.*Hx); HyH = 1./(Hy.*Hy);
C = ro.*0.5./(HxH+HyH); A = HxH.*C; B = HyH.*C; D = 1-ro;
A2 = 2.*A;
% формування крайових умов та початкового наближення :
e1 = clock();
U = zeros(Ny, Nx);
for i = 1 : Nx11               % D1 = {0<=x<0.3&0.3=y}
   U(Ny11,i) = mux11(X(i));
end
for i = Nx11 : Nx              % D3 = {0.3<=x<=1&0=y}
   U(1,i) = mux12(X(i));
end
for i = 1 : Nx21               % D5 = {0<=x<0.5&1=y}
   U(Ny,i) = mux21(X(i));
end
for j = 2 : Ny1                % D4 = {0.5<=x<=1&y=2-2*x}
   i = Nx-j+1; U(j,i) = muxy(X(i));
end

it  = 0; Fla = 0;
while ~Fla & it < itmax        % ітер.процес
   Fla = 1; it = it+1;

   for j = 2 : Ny11-1          % x=0.3&0<y<0.3
      i = Nx11; u = U(j,i);
      UU = A2.*(U(j,i+1)-Hx.*muy12(Y(j))) + ...
           B.*(U(j-1,i)+U(j+1,i)) + D.*u + C.*f(X(i),Y(j));
      if Fla
         Fla = abs(UU-u) < ep;
      end
      U(j,i) = UU;
      for i = Nx11+1 : Nx-j    % 0.3<x<1&0<y<0.3&y=2-2x 
         u = U(j,i);
         UU = A.*(U(j,i-1)+U(j,i+1)) + ...
              B.*(U(j-1,i)+U(j+1,i)) + D.*u + C.*f(X(i),Y(j));
         if Fla
            Fla = abs(UU-u) < ep;
         end
         U(j,i) = UU;
      end
   end

   j = Ny11;
   for i = Nx11+1 : Nx-j       % 0.3<x<1&y=0.3&y=2-2x
      u = U(j,i);
      UU = A.*(U(j,i-1)+U(j,i+1)) + ...
          B.*(U(j-1,i)+U(j+1,i)) + D.*u + C.*f(X(i),Y(j));
      if Fla
         Fla = abs(UU-u) < ep;
      end
      U(j,i) = UU;
   end
  
   for j = Ny11+1 : Ny1;       % 0=x&0.3<y<1&y=2-2x
      i = 1; u = U(j,i);
      UU = A2.*(U(j,i+1) - Hx.*muy11(Y(j))) + ...
           B.*(U(j-1,i)+U(j+1,i)) + D.*u + C.*f(X(i),Y(j));
      if Fla
         Fla = abs(UU-u) < ep;
      end
      U(j,i) = UU;
      for i = 2 : Nx-j         % 0<x<1&0.3<y<1&y=2-2x
         u = U(j,i);
         UU = A.*(U(j,i-1)+U(j,i+1)) + ...
              B.*(U(j-1,i)+U(j+1,i)) + D.*u + C.*f(X(i),Y(j));
         if Fla
            Fla = abs(UU-u) < ep;
         end
         U(j,i) = UU;
      end
   end

end

e2 = clock(); u = etime(e2,e1);
fprintf('Виконано ітерацій =%5u\n', it);
fprintf('Тривалість розрах.=%9.3f сек\n', u);

% побудова графіку розв'язку :
u = @(x,y)(15.*x.*(1-x).*y.*(1-y)); % точний розв'язок
UU = zeros(Ny,Nx);
for j = 1 : Ny
   for i = 1 : Nx
     if (j >= Ny11 | i >= Nx11) & i <= Nx-j+1
        UU(j,i) = u(X(i),Y(j));
     end
  end
end 
[Mx,My] = meshgrid(X,Y); clear X Y u
u = 5e-4; e1 = [-u,a+u,-u,b+u,-u,1+u];
surfl(Mx, My, U); axis(e1);
e2 = 'Різницевий розв’язок';
title(e2); xlabel('x'); ylabel('y'); zlabel('U');
shading interp; colormap(pink); view(40,20);
pause; close all;
contour(Mx, My, U, 15);
xlabel('x'); ylabel('y'); title(e2); colormap hsv;
pause; close all;

surfl(Mx, My, UU); axis(e1);
e2 = 'Точний розв’язок';
title(e2); xlabel('x'); ylabel('y'); zlabel('u');
shading interp; colormap(pink); view(40,20);
pause; close all;
contour(Mx, My, UU, 15);
xlabel('x'); ylabel('y'); title(e2); colormap hsv;
pause; close all;

UU = U - UU;
surfl(Mx, My, UU);
e2 = 'Похибка U-u';
title(e2); xlabel('x'); ylabel('y'); zlabel('U-u');
shading interp; colormap(pink); view(40,20);
e1 = norm(UU,inf); fprintf([e2,'  =%13.5e\n'], e1);
e1 = max(max(abs(UU))); fprintf('max(max|U-u|)=%13.5e\n', e1);
end

function ff = f(x, y)
% функція   f(x,y)
   ff = 30.*(x.*(1 - x) + y.*(1 - y));
end

function mu = mux11(x)
%  крайові умови на D1 = {(0<=x<=0.3)&(y=0.3)}
   mu = 3.15.*x.*(1 - x);
end

function mu = mux12(x)
%  крайові умови на D3 = {(0.3<=x<=1)&(y=0)}
   mu = 0;
end

function mu = mux21(x)
%  крайові умови на D5 = {(0<=x<=0.5)&(y=1)}
   mu = 0;
end

function mu = muxy(x)
% крайові умови на D4 = {(0.5<=x<=1)&(y=2-2x)}
   xx = 2.*x;
   mu = 15.*xx.*(1 - x).^2.*(xx - 1);
end

function mu = muy11(y)
%  крайові умови на D6 = {(0.3<y<1)&(x=0)}
   mu = 15.*y.*(1 - y);
end

function mu = muy12(y)
%  крайові умови на D2 = {(0<y<0.3)&(x=0.3)}
   mu = 6.*y.*(1 - y);
end