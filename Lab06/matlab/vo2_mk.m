function [X Y n] = vo2_mk(lambda, a, b, n, kfun, ffun)
% Розв'язок лiнiйного неоднорiдного iнтегрального рiвняння
% типу Вольтерри другого роду методом замiни iнтеграла
%  скiнченною сумою (квадратурнi формули методу трапецiй)
%                       x
%     y(x)=f(x)+lambda* S K(x,s)*y(s)ds,   x є [a,b].
%                       a
% вхiднi данi :
%  lambda - числовий параметр рiвняння;
%     a,b - інтервал визначення функції y(x);
%       n - число точок сiтки (n>1)
%    kfun - функцiя з описом ядра K(x,s);
%    ffun - функцiя з описом правої частини f(x);
% вихiднi данi :
%  X      - масив вузлiв сiтки;
%  Y      - масив наближеного розв'язку;
%  n      - розмірність масивiв X, Y.
if nargin<6
   error('VO2_MK:NotEnoughInputs',...
         'Not enough input arguments.  See VO2_MK.');
   return
end
if a >= b
   error('VO2_MK:error a >= b.','a >= b.  See VO2_MK.');
   return
end
if n < 2
  n = 2;
end
% формування вузлiв
X = linspace(a, b, n); h = X(2) - X(1);
Y = zeros(n, 1); P = Y;
%  Знаходження розв'язку (СЛАР з лівою трикутною матрицею):
%   Y(i) = (f(X(i))+lambda*SUM(k=1,i-1, P(k)*K(X(i),X(k))*Y(k) ))/
%            (1-lambda*P(i)*K(X(i),X(i))),  i=1,...,n.
for i = 1 : n
% формування коефiцiєнтiв P(j), j=1,...,i КФ
   P(1) = 0.5.*h;
   for j = 2 : i - 1
      P(j) = h;
   end
   P(i) = P(1);
   x = X(i); c1 = 0;
   for k = 1 : i - 1
      c1 = c1 + P(k).*kfun(x, X(k)).*Y(k);
   end
   c1 = ffun(x) + lambda.*c1;
   c2 = 1 - lambda.*P(i).*kfun(x, x);
   Y(i) = c1./c2;
end
return