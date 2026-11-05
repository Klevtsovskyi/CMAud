function U =  fdm_P10(p, f, u0, Ga, Gb, v, X, J, tau)
% Наближений розв'язок крайової задачі :
% du(x,t)/dt = d^2u(x,t)/dx^2 - p(x,t)*u + f(x,t),
%            p(x)>=0, x є (a,b),
%         u(x,0) = u0(x),  x є [a,b],
% КУ I роду      va=0 :                   u(a,t) = мю(t)
% КУ II-III роду va=1 : -du(a,t)/dx + g(t)u(a,t) = мю(t)
% КУ I роду      vb=0 :                   u(b,t) = ню(t)
% КУ II-III роду vb=1 :  du(b,t)/dx + q(t)u(b,t) = ню(t)
%        (явна схема  sigm = 0)
% Вхiднi данi :
%      p - функцiя з описом p(x,t)
%      f - функцiя з описом f(x,t)
%     u0 - функцiя з описом u0(x)
%     Ga - вектор-функція [   1, мю(t)] при va=0
%                         [g(t), мю(t)] при va=1
%     Gb - вектор-функція [   1, ню(t)] при vb=0
%                         [q(t), ню(t)] при vb=1
%      v - [va, vb] - масив значень виду КУ в x=a,b
%      X - масив рівномірної сітки по x a=X(1),b=X(end)
%      J - масив з номерами шарів для запам'ятовування розв'язку
%          max(J(:)) = M - гранична точка сітки по t, T=M*tau
%    tau - крок сітки по змінній t
% Вихiднi данi :
%     U  - масив значень наближеного розв'язку
%           U(j,k) = u(X(k),tau*J(j))

% параметри різницевої схеми :
n = length(X); nm = n-1; h = X(2)-X(1); hh = h.*h;
M = length(J); Y = X; U = zeros(M,n); M = max(J); 
if 2*tau >= hh
   warning('FDM_P10: не виконується умова стійкості!');   
end
% обчислення додаткових констант
tau05 = 0.5.*tau;
gam = tau./hh; gam2 = 2.*gam; gamt = 2.*tau./h; g12 = 1-gam2;
% формування початкових умов при t = 0 :
for k = 1 : n
   Y(k) = u0(X(k));
end
U(1,:) = Y; nJ = 2; % запам'ятовуємо для побудови 3D-графіку
for j = 1 : M
   t = j.*tau; tt = t - tau05;
% знаходження наближеного розв'язку Y(:) на шарі j:
   yc = Y(1);
   if v(1) == 0
      gr = Ga(t); Y(1) = gm(2);
   else
      gr = Ga(tt);
      x = X(1); pt = tau05.*p(x,tt); P = 1 + pt;
      F = tau.*f(x,tt) + gamt.*gr(2);
      x = g12 - pt - gamt.*gr(1);
      Y(1) = (x.*yc + gam2.*Y(2) + F)./P;
   end
   for k = 2 : nm
      x = X(k); yl = yc; yc = Y(k);
      P = tau.*p(x,tt); F = tau.*f(x,tt);
      Y(k) = gam.*(yl+Y(k+1)) + (g12-P).*yc + F;
   end
   if v(2) == 0
      qr = Gb(t); Y(n) = qn(2);
   else
      qr = Gb(tt);
      x = X(n); pt = tau05.*p(x,tt); P = 1 + pt;
      F = tau.*f(x,tt) + gamt.*qr(2);
      x = g12 - pt - gamt.*qr(1);
      Y(n) = (gam2.*yc + x.*Y(n) + F)./P;
   end
   if j == J(nJ)
      U(nJ,:) = Y; nJ = nJ+1; % для побудови 3D-графіку
   end
end
end