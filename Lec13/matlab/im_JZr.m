function [U, it] = im_JZr(m, f, muh, muv, X, Y, ro, ep, itmax)
% m=0 м.Якобі, m=1 м.Зейделя, m=2 м.релаксації
   % параметри різницевої схеми :
   Nx = length(X); Nx1 = Nx - 1; Hx = X(2) - X(1);
   Ny = length(Y); Ny1 = Ny - 1; Hy = Y(2) - Y(1);
   % обчислення констант розрахункових формул :
   if m == 2
      Hx2 = 1./(Hx.*Hx); Hy2 = 1./(Hy.*Hy); H = 2.*(Hx2+Hy2); 
      C = ro./H; A = Hx2.*C; B = Hy2.*C; D = 1-ro;
   else
      Hx2 = Hx.*Hx; Hy2 = Hy.*Hy; H = 2.*(Hx2+Hy2); 
      A = Hy2./H; B = Hx2./H; C = Hx2.*Hy2./H; D = 0;
   end
   % формування граничних умов і початкового наближення :
   U = zeros(Ny, Nx);
   for i = 1 : Nx         % КУ на "горизонталях" y=0, y=b
      [c1, c2] = muh(X(i)); U(1,i) = c1; U(Ny,i) = c2;
   end
   for j = 2 : Ny1        % КУ на "вертикалях" x=0, x=a  
       [c1, c2] = muv(Y(j));
       % + лінійна інтерполяція по горизонталі
       U(j,:) = linspace(c1, c2, Nx);
   end
   clear c1 c2

   it  = 0; Fla = 0;
   while ~Fla & it < itmax        % ітер.процес
      Fla = 1; it = it+1;
      if m == 0
         Uj = U(1,:);          
      end      
      for j = 2 : Ny1             % рух "знизу вгору"
         jm = j-1; y = Y(j);
         if m == 0
            Ujm = Uj; Uj = U(j,:);
         end
         for i = 2 : Nx1          % рух "зліва направо"
            if m == 0
               old0 = Ujm(i);  oldm = Uj(i-1);
            else
               old0 = U(jm,i); oldm = U(j,i-1);
            end
            old = U(j,i);
            new = A.*(oldm + U(j,i+1)) + ...
                  B.*(old0 + U(j+1,i)) + C.*f(X(i),y);
            if m == 2
                new = new + D.*old;
            end
            if Fla
               Fla = abs(new-old) < ep;
            end
            U(j,i) = new;
         end
      end
   end
end