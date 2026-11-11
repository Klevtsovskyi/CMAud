function U =  fdm_P1n2(fi, fdf, u0, mn, X, J, tau, ep, smax)
% Наближений розв'язок крайової задачі :
% d(фі(u))/dt = u"(x,t) + f(u),
%           x є (a,b),
%      u(x,0) = u0(x),  x є [a,b],
%    u(a,t) = мю(t),  u(b,t) = ню(t).
% Вхiднi данi :
%     fi - функцiя з описом фі(u)
%    fdf - функцiя з описом [f(u), d(f(u))/du, d(фі(u))/du]
%     u0 - функцiя з описом u0(x)
%     mn - вектор-функція з описом [мю(t), ню(t)]
%      X - масив рівномірної сітки по x a=X(1),b=X(end)
%      J - масив з номерами шарів для запам'ятовування розв'язку
%          max(J(:)) = M - гранична точка сітки по t, T=M*tau
%    tau - крок сітки по змінній t
%     ep - точність (0<ep<1)
%   smax - max кількість ітерацій
% Вихiднi данi :
%     U  - масив значень наближеного розв'язку
%           U(j,k) = u(X(k),tau*J(j))

% параметри різницевої схеми :
n = length(X); nm = n-1; h = X(2)-X(1); hh = h.*h;
M = length(J); Y = X; U = zeros(M,n); M = max(J); 
% обчислення додаткових констант
gam = tau./hh; gam2 = 2.*gam;
% формування початкових умов при t = 0 :
for i = 1 : n
   Y(i) = u0(X(i));
end
U(1,:) = Y; nJ = 2; % запам'ятовуємо для побудови 3D-графіку
% A,C,B,F - елементи СЛАР
A = ones(1,n); F = A; C = A;
A(1) = 0; A(n) = 0; A = gam.*A; B = A;
for j = 1 : M
   t = j.*tau; Yp = Y; m_n = mn(t);  
   for s = 1 : smax
      Ys = Y;
      % формуємо решту коефіцієнтів СЛАР на шарі j:
      F(1) = m_n(1) - Ys(1); F(n) = m_n(2) - Ys(n);
      for i = 2 : nm
         x = X(i); y = Ys(i); f_df = fdf(y);
         C(i) = gam2 - tau.*f_df(2) + f_df(3);
         F(i) = gam.*(Ys(i-1)-2.*y+Ys(i+1)) + ...
                tau.*f_df(1) + fi(Yp(i)) - fi(y);
      end
      % знаходження Y на s-й ітерації для шару j:
      [Y AL] = m_progn(A, C, B, F);
      x = max(abs(Y)); Y = Ys + Y;
      if x <= ep  % закінчуємо ітерацію по s
         break;
      end
   end
   if j == J(nJ)
      U(nJ,:) = Y; nJ = nJ+1; % для побудови 3D-графіку
   end
end
end