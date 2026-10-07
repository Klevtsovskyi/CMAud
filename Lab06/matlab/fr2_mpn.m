function [X Y n] = fr2_mpn(lambda, a, b, nkf, n, e, mki, kfun, ffun)
% Розв'язок лiнiйного неоднорiдного iнтегрального рiвняння
% типу Фредгольма другого роду методом послiдовних наближень
%                       b
%     y(x)=f(x)+lambda* S K(x,s)*y(s)ds
%                       a
% вхiднi данi :
%  lambda  - числовий параметр рiвняння;
%      a   - нижня границя iнтеграла;
%      b   - верхня границя iнтеграла;
%      nkf - 0 = КФ методу трапецiй,
%            1 = КФ методу Сiмпсона;
%       n  - число точок сiтки. Повинно задовольняти умовi:
%            n > 1,                 при nkf=0,
%            n > 2 & mod(n,2) == 1, при nkf=1;
%       e  - точнiсть (e>0);
%     mki  - максимальна кiлькiсть крокiв iтерацiй;
%     kfun - функцiя з описом ядра K(x,s);
%     ffun - функцiя з описом правої частини f(x);
% вихiднi данi:
%   X      - масив вузлiв сiтки;
%   Y      - масив наближеного розв'язку;
%   n      - розмірність масивiв X, Y.

if nargin<9
   error('FR2_MPN:NotEnoughInputs',...
         'Not enough input arguments.  See FR2_MPN.');
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
            n = n + 1;
         end
      end
   otherwise
      error('FR2_MPN:Unknown method.','Unknown method.  See FR2_MPN.');
      return
end
if a >= b
   error('FR2_MPN:error a >= b.','a >= b.  See FR2_MPN.');
   return
end
if e <= 0
   error('FR2_MPN:error e <= 0.','e <= 0.  See FR2_MPN.');
   return
end
% формування вузлiв
X = linspace(a, b, n); h = X(2) - X(1); n1 = n - 1;
% формування коефiцiєнтiв P(i), i=1,...,n  КФ
P = zeros(n, 1); FI = P; Y = FI; AK = zeros(n);
switch nkf
   case 0
      P(1) = 0.5 .* h;
      for i=2:n1
         P(i) = h;
      end
      P(n) = P(1);
   case 1
      P(1) = h ./ 3; c2 = 2 .* P(1);  c1 = 2 .* c2;
      for i = 2 : 2 : n1
         P(i) = c1;  P(i+1) = c2;
      end
      P(n) = P(1);
end
% формування матрицi AK(i,j) i вектора FI(i),
% де AK(i,j)=lambda*P(j)*K(X(i),X(j)),
%    FI(i)  = f(X(i)),  i=1,...,n; j=1,...,n
for i = 1 : n
   x = X(i); ff = ffun(x); FI(i) = ff;  Y(i) = ff;
   for j = 1 : n
      AK(i, j) = lambda .* P(j) .* kfun(x,X(j));
   end
end
% послiдовні наближення
k = 0; F = 1;
while F & k < mki
   k = k + 1; F = 0; FIi = AK * FI;
   for i = 1 : n
      ff  = FIi(i); Y(i) = Y(i) + ff; FI(i) = ff;
      if F == 0
         F = abs(ff) > e;
      end
   end
end
return