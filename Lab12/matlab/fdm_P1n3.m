function U =  fdm_P1n3(K, f, u0, mn, X, J, tau)
% Наближений розв'язок крайової задачі :
% du(x,t)/dt = d(K(u)*u')/dx + f(u),
%       K > 0, x є (a,b),
%     u(x,0) = u0(x),  x є [a,b],
%   u(a,t) = мю(t), u(b,t) = ню(t).
% Вхiднi данi :
%      K - функцiя з описом K(u)
%      f - функцiя з описом f(u)
%     u0 - функцiя з описом u0(x)
%     mn - вектор-функція з описом [мю(t), ню(t)]
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
% обчислення додаткових констант
gam = tau./hh; gam05 = 0.5.*gam;
tau05 = 0.5.*tau; tau2 = 2.*tau;
% формування початкових умов при t = 0 :
for i = 1 : n
   Y(i) = u0(X(i));
end
U(1,:) = Y; nJ = 2; % запам'ятовуємо для побудови 3D-графіку
% A,C,B,F - елементи СЛАР
A = zeros(1,n); B = A; F = A; C = ones(1,n);
for j = 1 : M
   t = (j-0.5).*tau;
   % формуємо СЛАР на шарі j+0.5:
   m_n = mn(t); F(1) = m_n(1); F(n) = m_n(2);  
   ki = K(Y(1)); kip = K(Y(2));
   for i = 2 : nm
      ip = i+1; y = Y(i);
      kim = ki; ki = kip; kip = K(Y(ip));
      ai = gam05.*(kim + ki); A(i) = ai;
      bi = gam05.*(ki + kip); B(i) = bi;
      C(i) = ai + bi + 2;
      F(i) = 2.*y + tau.*f(y);
   end
   % знаходження Y для шару j+0.5:
   [Yt, AL, er] = m_progm(A, C, B, F);
   % формуємо СЛАР на шарі j+1:
   t = t + tau05;
   m_n = mn(t); F(1) = m_n(1); F(n) = m_n(2);  
   ki = K(Yt(1)); kip = K(Yt(2));
   for i = 2 : nm
      ip = i+1; y = Yt(i);
      kim = ki; ki = kip; kip = K(Yt(ip));
      ai = gam05.*(kim + ki); A(i) = ai;
      bi = gam05.*(ki + kip); B(i) = bi;
      x = ai + bi; C(i) = 2 + x; x = 2 - x;
      F(i) = ai.*Y(i-1)+x.*Y(i)+bi.*Y(ip) + tau2.*f(y);
   end
   % знаходження Y для шару j+1:
   [Y, AL, er] = m_progm(A, C, B, F);
   if j == J(nJ)
      U(nJ,:) = Y; nJ = nJ+1; % для побудови 3D-графіку
   end
end
end