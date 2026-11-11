function view_solu(X, J, Y, U, ti1, z1, ti2, z2, txtp, ti3, z3, s)
close all;
n = length(X); m = length(J);
mU = inf; MU = -inf;
for i = 1 : n
   x = X(i);
   for j = 1 : m
      de = U(j,i);
      if de > MU
         MU = de;
      elseif de < mU
         mU = de;
      end
   end
end
[Mx, Mt] = meshgrid(X, J); % побудова сіток для 3D-графіків
% виведення графіка розв'язку Y
de = [X(1), X(end), J(1), J(m), mU-0.1, MU+0.1];
figure(1); mesh(Mx, Mt, Y); axis(de);
shading interp, colormap cool;
xlabel('x'), ylabel('t'), zlabel(z1);
title(ti1);
if s, pause(s);
else, pause;
end
% виведення графіка U
figure(2); mesh(Mx, Mt, U); axis(de);
shading interp, colormap cool;
xlabel('x'); ylabel('t'); zlabel(z2);
title(ti2);
if s, pause(s);
else, pause
end
% виведення графіка Y - U
YU = Y - U;
C = norm(YU,inf);
fprintf('%s ||Y-U||c=%13.5e\n',txtp,C);
C = max(max(abs(YU)));
fprintf('max(max(|Y-U|)) =%13.5e\n',C);
figure(3); mesh(Mx, Mt, YU); view(-145,35);
shading interp, colormap jet;
xlabel('x'); ylabel('t'); zlabel(z3);
title(ti3);
end
