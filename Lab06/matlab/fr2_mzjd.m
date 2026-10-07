function [X Y n] = fr2_mzjd(lambda,a,b,nkf,n,kfun,ffun)
% Розв'язок лiнiйного неоднорiдного iнтегрального рiвняння
% типу Фредгольма другого роду методом замiни ядра виродженим
%                       b
%     y(x)=f(x)+lambda* S K(x,s)*y(s)ds,  
%                       a
%      де K(x,s)=SUM(k=1,m,Ak(x)*Bk(s))
% вхiднi данi :
%  lambda - числовий параметр рiвняння;
%       a - нижня границя iнтеграла;
%       b - верхня границя iнтеграла;
%     nkf - 0 = КФ методу трапецiй,
%           1 = КФ методу Сiмпсона 
%            (для чисельного обчислення iнтегралiв);
%       n - число точок сiтки. Повинно задовольняти умовi:
%           n > 1,                     при nkf=0,
%           n > 2 & mod(n,2) == 1,     при nkf=1;
%    kfun - функцiя з описом ядра K(x,s);
%    ffun - функцiя з описом правої частини f(x);
% вихiднi данi :
%  X      - масив вузлiв сiтки;
%  Y      - масив наближеного розв'язку;
%  n      - розмірність масивiв X, Y.
global m
if nargin<7
   error('FR2_MZJD:NotEnoughInputs',...
         'Not enough input arguments.  See FR2_MZJD.');
   return
end
if a >= b
   error('FR2_MZJD:error a >= b.','a >= b.  See FR2_MZJD.');
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
      error('FR2_MZJD:Unknown method.',...
            'Unknown method.  See FR2_MZJD.');
      return
end
% формування вузлiв
X = linspace(a, b, n); h = X(2) - X(1); n1 = n - 1;
% формування коефiцiєнтiв P(i), i=1,...,n КФ
P = zeros(n,1); FF = P; Y = P;
switch nkf
   case 0
      P(1) = 0.5.*h;
      for i = 2 : n1
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
% формування допомiжних матриць AK, BK i вектора FF,
%    де AK(k,i)=Ak(X(i)), BK(k,i)=Bk(X(i)),
%         FF(i)=f(X(i)),  i=1,...,n.
for i = 1 : n
  x = X(i);
  [AK(:, i), BK(:, i)] = kfun(x); FF(i) = ffun(x);
end
% формування матрицi A i вектора F  СЛАР  A*C=F,
%    де  A(i,j)=if(i=j,1,0)-lambda*SUM(k=1,n, P(k)*AK(j,k)*BK(i,k) ),
%        F(i)  =SUM(k=1,n, P(k)*FF(k)*BK(i,k) ),  i=1,...,n
F = zeros(m,1); A = zeros(m,m);
for i = 1 : m
  for j = 1 : m
     aa = 0;
     for k = 1 : n    
        aa = aa + P(k).*AK(j, k).*BK(i, k);
     end
     aa = -lambda.*aa;
     if i  ==  j
        aa = 1 + aa;
     end
     A(i,j) = aa;
  end
  aa = 0;
  for k = 1 : n
    aa = aa + P(k).*FF(k).*BK(i, k);
  end
  F(i) = aa;
end
% розв'язок СЛАР    A*C=F
C = A\F;
% розв'язок у виглядi
% y(x)=f(x)+lambda*SUM(k=1,m, C(k)*Ak(x) )
for i = 1 : n
   Y(i) = FF(i) + lambda.*dot(C, AK(:, i));
end
return