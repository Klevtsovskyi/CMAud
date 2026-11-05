function U =  fdm_P1s(p, f, u0, Ga, Gb, v, X, J, tau, sigm)
% Наближений розв'язок крайової задачі :
% du(x,t)/dt = d^2u(x,t)/dx^2 - p(x,t)*u + f(x,t),
%            p(x)>=0, x є (a,b),
%         u(x,0) = u0(x),  x є [a,b],
% КУ I роду      va=0 :                   u(a,t) = мю(t)
% КУ II-III роду va=1 : -du(a,t)/dx + g(t)u(a,t) = мю(t)
% КУ I роду      vb=0 :                   u(b,t) = ню(t)
% КУ II-III роду vb=1 :  du(b,t)/dx + q(t)u(b,t) = ню(t)
%        (неявна схема  0 <sigm<= 1)
% Вхiднi данi :
%      p - функцiя з описом p(x,t)
%      f - функцiя з описом f(x,t)
%     u0 - функцiя з описом u0(x)
%     Ga - вектор-функція [   1, мю(t)] при va=0
%                         [g(t), мю(t)] при va=1
%     Gb - вектор-функція [   1, ню(t)] при vb=0
%                         [q(t), ню(t)] при vb=1
%      v - [va, vb] - масив означень виду КУ в x=a,b
%      X - масив рівномірної сітки по x a=X(1),b=X(end)
%      J - масив з номерами шарів для запам'ятовування розв'язку
%          max(J(:)) = M - гранична точка сітки по t, T=M*tau
%    tau - крок сітки по змінній t
%   sigm - вага РС (0<sigm<=1)
% Вихiднi данi :
%     U  - масив значень наближеного розв'язку
%           U(j,k) = u(X(k),tau*J(j))

% параметри різницевої схеми :
n = length(X); nm = n-1; h = X(2)-X(1); hh = h.*h;
M = length(J); Y = X; U = zeros(M,n); M = max(J); 
if sigm <= 0 | sigm > 1
   sigm = 0.5;
   warning('FDM_P1S: значення ваги sigm змінено на 0.5!');   
end
% обчислення додаткових констант
tau05 = 0.5.*tau; s1 = 1-sigm;
gam = tau./hh; gam2 = 2.*gam; gamt = 2.*tau./h;
gams = s1.*gam; gam2s =s1.*gam2;
gam = sigm.*gam; gam2 =sigm.*gam2;
g12p = 1+gam2; g12m = 1-gam2s;
sgt = sigm.*gamt; s1gt = s1.*gamt;
st = sigm.*tau; s1t = s1.*tau;
% формування початкових умов при t = 0 :
for k = 1 : n
   Y(k) = u0(X(k));
end
U(1,:) = Y; nJ = 2; % запам'ятовуємо для побудови 3D-графіку
% формуємо A,C,B,F - елементи СЛАР
A = ones(1,n); C = A; F = A; A = gam.*A; B = A; 
A(1) = 0;
if v(1) == 0
   B(1) = 0;
else
   B(1) = gam2;    
end
if v(2) == 0
   A(n) = 0;
else
   A(n) = gam2;    
end
B(n) = 0;

for j = 1 : M
   t = j.*tau; tt = t - tau05;
% формуємо коефіцієнти СЛАР для шару t=t(j):
   if v(1) == 0
      gr = Ga(t); F(1) = gr(2);
   else
      gr = Ga(tt); x = X(1); pt = tau05.*p(x,tt);
      C(1) = g12p+pt+sgt.*gr(1);
      F(1) = (g12m-pt-s1gt.*gr(1)).*Y(1) + gam2s.*Y(2) ...
             + tau.*f(x,tt) + gamt.*gr(2);
   end
   for k = 2 : nm
      x = X(k); pt = p(x,tt);
      C(k) = g12p + st.*pt;
      F(k) = gams.*(Y(k-1)+Y(k+1)) + (g12m-s1t.*pt).*Y(k) ...
             + tau.*f(x,tt);
   end
   if v(2) == 0
      qr = Gb(t); F(n) = qr(2);
   else
      qr = Gb(tt); x = X(n); pt = tau05.*p(x,tt);
      C(n) = g12p + pt + sgt.*qr(1);
      F(n) = gam2s.*Y(nm) + (g12m-pt-s1gt.*qr(1)).*Y(n) ...
             + tau.*f(x,tt) + gamt.*qr(2);
   end
% знаходження наближеного розв'язку Y(:) на шарі t=t(j):
   [Y AL, er] = m_progm(A, C, B, F);
   if j == J(nJ)
      U(nJ,:) = Y; nJ = nJ+1; % для побудови 3D-графіку
   end
end
end