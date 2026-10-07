function  [X Y n]  =  fr2_mk(lambda, a, b, nkf, n, kfun, ffun)
% Розв'язок лiнiйного неоднорiдного iнтегрального рiвняння
% типу Фредгольма другого роду методом замiни iнтеграла
%  скiнченною сумою
%                       b
%     y(x)=f(x)+lambda* S K(x,s)*y(s)ds
%                       a
% вхiднi данi :
%  lambda - числовий параметр рiвняння;
%  a      - нижня границя iнтеграла;
%  b      - верхня границя iнтеграла;
%  nkf    - 0 = КФ методу трапецiй,
%           1 = КФ методу Сiмпсона;
%  n      - число точок сiтки. Повинно задовольняти умовi:
%           n > 1,                  при nkf=0,
%           n > 2 & mod(n,2) == 1,  при nkf=1;
%  kfun   - функцiя з описом ядра K(x,s);
%  ffun   - функцiя з описом правої частини f(x);
% вихiднi данi:
%  X      - масив вузлiв сiтки;
%  Y      - масив наближеного розв'язку;
%  n      - розмірність масивiв X, Y.

if nargin<7
   error('FR2_MK:NotEnoughInputs',...
         'Not enough input arguments.  See FR2_MK.');
   return
end
if a  >=  b
   error('FR2_MK:error a >= b.',...
          'a >= b.  See FR2_MK.');
   return
end
switch nkf
   case 0
      if n < 2
         n = 2;
      end
   case 1
      if n < 3
         n = 3;
      else
         if mod(n,2) == 0
            n = n+1;
         end
      end
   otherwise
      error('FR2_MK:Unknown method.',...
            'Unknown method.  See FR2_MK.');
      return
end
% формування вузлiв
X = linspace(a,b,n); h = X(2) - X(1); n1 = n-1;
% формування коефiцiєнтiв P(i), i=1,...,n КФ
P = zeros(n,1);
switch nkf
   case 0
      P(1) = 0.5.*h;
      for i=2:n1
         P(i) = h;
      end
      P(n) = P(1);
   case 1
      P(1) = h./3; c2 = 2.*P(1);  c1 = 2.*c2;
      for i=2:2:n1
         P(i) = c1; P(i+1) = c2;
      end
      P(n) = P(1);
end
% формування матрицi A(i,k) i вектора правої частини F(i) для СЛАР
%   Y(i)-lambda*SUM(k=1,n, P(k)*K(X(i),X(k))*Y(k) ) = f(X(i)),
%   i=1,...,n
A = zeros(n, n);  F = P;
for i = 1 : n
   x = X(i);  F(i) = ffun(x);
   for k = 1 : n
      A(i, k) = -lambda .* P(k) .* kfun(x, X(k));
   end
   A(i, i) = 1 + A(i, i);
end
% розв'язок СЛАР    A*Y=F
Y = A\F;
return