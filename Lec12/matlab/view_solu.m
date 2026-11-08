function view_solu(X1, X2, Y, z1, ti1, pnt, s, U, z2, ti2, txtp, z3, ti3)
N1 = length(X1); N2 = length(X2);
mY = inf; MY = -inf;  % mY=min(Y), MY=max(Y)
for k2 = 1 : N2
   for k1 = 1 : N1
      de = Y(k2,k1);
      if de > MY
         MY = de;
      elseif de < mY
         mY = de;
      end
   end
end

if pnt   %  при необхідності - друк
   for k2 = N2:-1:1
      fprintf('%11.2e',Y(k2,:));
      fprintf('\n');
   end
end

[Mx1, Mx2] = meshgrid(X1, X2); % побудова сіток для 3D-графіків
% виведення графіка розв'язку Y
de = [X1(1), X1(N1), X2(1), X2(N2), mY-0.1, MY+0.1];
mesh(Mx1, Mx2, Y); axis(de);
shading interp, colormap cool;
xlabel('x1'), ylabel('x2'), zlabel(z1); title(ti1);
if s, pause(s);
else, pause;
end
close all

if nargin == 7
   return
end

% виведення графіка U
mesh(Mx1, Mx2, U); axis(de);
shading interp, colormap cool;
xlabel('x1'); ylabel('x2'); zlabel(z2); title(ti2);
if s, pause(s);
else, pause
end
close all
% виведення графіка Y - U
YU = Y - U;
C = norm(YU,inf); 
fprintf('%s ||Y-U||c=%13.5e\n',txtp,C);
C = max(max(abs(YU))); 
fprintf('max(max(|Y-U|)) =%13.5e\n',C);
mesh(Mx1, Mx2, YU); %view(-145,35);
shading interp, colormap jet;
xlabel('x1'); ylabel('x2'); zlabel(z3); title(ti3);
end