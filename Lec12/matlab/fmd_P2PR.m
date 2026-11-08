function Y = fmd_P2PR(f, muh, muv, u0, X1, X2, tau, M, J)
% Розв'язування лінійної мішаної задачі для двовимірного
% рівняння параболічного типу з КУ І роду за допомогою ППС.
% Вхiднi данi :
%      f - функцiя з описом f(x1,x2,t)
%    muh - функцiя мю на "горизонталі" [мю(x1,X2(1),t),мю(x1,X2(N2),t)]
%    muv - функцiя мю на "вертикалі"   [мю(X1(1),x2,t),мю(X1(N1),x2,t)]
%     u0 - функцiя з описом u0(x1,x2)
%  X1,X2 - масиви рівномірних сіток по змінним x1,x2
%    tau - крок сітки по змінній t
%      M - гранична точка сітки по t (T=M*tau)
%      J - масив з номерами шарів для запиcу розв'язку на МД.
%          Має бути J(1)<J(2)<...<J(end) або пустий. 
% Вихiднi данi :
%     Y  - масив значень наближеного розв'язку
%           Y(k2,k1) = u(X1(:),X2(:)) для j = M

% параметри різницевої схеми :
N1 = length(X1); N1m = N1 - 1; h1 = X1(2) - X1(1);
N2 = length(X2); N2m = N2 - 1; h2 = X2(2) - X2(1);
% обчислення додаткових констант :
tau05 = 0.5.*tau;
h12 = h1.*h1; ga1 = tau05./h12; s1 = 2.*ga1;
c1 = 1 + s1; s1 = 1 - s1;
h22 = h2.*h2; ga2 = tau05./h22; s2 = 2.*ga2;
c2 = 1 + s2; s2 = 1 - s2;
% формування елементів матриці СЛАР незалежних від j :
C1 = ones(1,N1); F1 = C1; A1 = [0, ga1.*C1(2:N1m), 0];
B1 = A1; C1(2:N1m) = c1;
C2 = ones(1,N2); F2 = C2; A2 = [0, ga2.*C2(2:N2m), 0];
B2 = A2; C2(2:N2m) = c2;
% формування початкових умов при t=0 :
Y = zeros(N2,N1); Y1 = Y;
for k2 = 1 : N2
   x2 = X2(k2);
   for k1 = 1 : N1
      Y(k2,k1) = u0(X1(k1),x2);
   end
end
nempJ = ~isempty(J);
if nempJ
   jout = 1; name = 'P2PR_';
   if J(jout) == 0
      namef = sprintf('%s%u.mat', name, J(jout));
      save(namef, 'Y'); jout = jout + 1;
   end
end
for j = 1 : M
% знаходження наближеного розв'язку Y на шарі j-1/2 :
   t = j.*tau; tt = t - tau05; tm = (j-1).*tau; 
% обчислення крайових умов в точках X2(1),X2(N2) :
   for k1 = 1 : N1
      [c1, c2] = muh(X1(k1),tt); Y1(1,k1) = c1; Y1(N2,k1) = c2;
   end
   for k2 = 2 : N2m
      x2 = X2(k2); x2m = X2(k2-1); x2p = X2(k2+1);
% формування F1 :
      [c1_m,c2_m] = muv(x2m,tm); [c1m,c2m] = muv(x2m,t);
      [c1_, c2_]  = muv(x2,tm);  [c1, c2]  = muv(x2,t);
      [c1_p,c2_p] = muv(x2p,tm); [c1p,c2p] = muv(x2p,t);
      F1(1) = 0.5.*((c1_+c1) + ga2.* ...
              ((c1_m-c1m)-2.*(c1_-c1)+(c1_p-c1p)));
      F1(N1)= 0.5.*((c2_+c2) + ga2.* ...
              ((c2_m-c2m)-2.*(c2_-c2)+(c2_p-c2p)));;
      for k1 = 2 : N1m
         F1(k1) = ga2.*(Y(k2,k1-1) + Y(k2,k1+1)) + ...
                  s2.*Y(k2,k1) + tau05.*f(X1(k1),x2,tt);
      end
% розв'язок СЛАР, знаходження Y(k2,:) :
      [y, alfa, err]  =  m_progm(A1, C1, B1, F1);
      Y1(k2,:) = y;
   end
   clear y alfa
% знаходження наближеного розв'язку Y на шарі j :
% обчислення крайових умов в точках X1(1),X1(N1) :
   for k2 = 1 : N2
      [c1, c2] = muv(X2(k2),t); Y(k2,1) = c1; Y(k2,N1) = c2;
   end
   for k1 = 2 : N1m
      x1 = X1(k1); k1m = k1 - 1; k1p = k1 + 1;
% формування F2 :
      [c1, c2] = muh(x1,t); F2(1) = c1; F2(N2) = c2;
      for k2 = 2 : N2m
         F2(k2) = ga1.*(Y1(k2,k1m) + Y1(k2,k1p)) + ...
                  s1.*Y1(k2,k1) + tau05.*f(x1,X2(k2),tt);
      end
% розв'язок СЛАР, знаходження Y(:,k1) :
      [y, alfa, err]  =  m_progm(A2, C2, B2, F2);
      Y(:,k1) = y';
   end
   clear y alfa
   if nempJ
      if J(jout) == j
         namef = sprintf('%s%u.mat', name, J(jout));
         save(namef, 'Y'); jout = jout + 1;
      end
   end
end
end