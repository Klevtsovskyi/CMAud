function U =  fdm_P1n1(K, cf, u0, mn, X, J, tau, ep, smax)
% Наближений розв'язок крайової задачі :
% c(x,t,u)*du(x,t)/dt = d(K(x,t,u)*u')/dx + f(x,t,u),
%            c, K > 0, x є (a,b),
%         u(x,0) = u0(x),  x є [a,b],
%    u(a,t) = мю(t),  u(b,t) = ню(t).
% Вхiднi данi :
%      K - функцiя з описом K(x,t,u)
%     cf - вектор-функцiя з описом [c(x,t,u), f(x,t,u)]
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
gam = tau./hh; gam05 = 0.5.*gam;
% формування початкових умов при t = 0 :
for i = 1 : n
   Y(i) = u0(X(i));
end
U(1,:) = Y; nJ = 2; % запам'ятовуємо для побудови 3D-графіку
% A,C,B,F - елементи СЛАР
A = zeros(1,n); B = A; F = A; C = ones(1,n);  
for j = 1 : M
   t = j.*tau; Yp = Y; 
   % формуємо F(1),F(n) для СЛАР на шарі j:
   m_n = mn(t); F(1) = m_n(1); F(n) = m_n(2);  
   for s = 1 : smax
      Ys = Y;
      % формуємо решту коефіцієнтів СЛАР на шарі j:
      ki = K(X(1),t,Y(1)); kip = K(X(2),t,Y(2));
      for i = 2 : nm
         ip = i+1; x = X(i); y = Y(i); c_f = cf(x,t,y);
         kim = ki; ki = kip; kip = K(X(ip),t,Y(ip));
         ai = gam05.*(kim + ki); A(i) = ai;
         bi = gam05.*(ki + kip); B(i) = bi;
         C(i) = ai + bi + c_f(1);
         F(i) = c_f(1).*Yp(i) + tau.*c_f(2);
      end
      % знаходження Y на s-й ітерації для шару j:
      [Y AL, er] = m_progm(A, C, B, F);
      if max(abs(Y-Ys)) <= ep  % закінчуємо ітерацію по s
         break;
      end
   end
   if j == J(nJ)
      U(nJ,:) = Y; nJ = nJ+1; % для побудови 3D-графіку
   end
end
end