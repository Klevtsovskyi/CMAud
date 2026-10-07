function [x Y] = fr2_mcol(lambda, X, kfun, ffun)
% Розв'язок лiнiйного неоднорiдного iнтегрального рiвняння
% типу Фредгольма другого роду методом колокацій
%                       b
%     y(x)=f(x)+lambda* S K(x,s)*y(s)ds
%                       a
% вхiднi данi :
%  lambda - числовий параметр рiвняння;
%  X      - сплайнова сітка a=X(1)<X(2)<...<X(end-1)<X(end)=b,
%            де [a,b] - інтервал інтегрування;
%  kfun   - функцiя з описом ядра K(x,s);
%  ffun   - функцiя з описом правої частини f(x);
% вихiднi данi :
%  x      - масив колокаційної сiтки x(i)=0.5*(X(i)+X(i+1)), i=1,...,n-1;
%  Y      - масив наближеного розв'язку в точках колокації;

if nargin<4
   error('FR2_MCOL:NotEnoughInputs',...
         'Not enough input arguments.  See FR2_MCOL.');
   return
end
n = length(X);
if n < 2
   n = 2;
end
% формування колокаційної сiтки
N = n-1; x = zeros(1, N); h = x; I = 1 : N;
for i = I
   j = i+1; 
   x(i) = 0.5.*(X(i)+X(j));
   h(i) = X(j)-X(i);
end
% формування матрицi A(i,j) i вектора правої частини F(i) для СЛАР
%   Y(i)-lambda*SUM(j=1,N,(X(j+1)-X(j))*K(x(i),x(j))*Y(j)) = f(x(i)),
%   i=1,...,N
A = zeros(N, N);  F = zeros(N, 1);
for i = I
   xi = x(i);  F(i) = ffun(xi);
   for j = I
      A(i, j) = -lambda.*h(j).*kfun(xi, x(j));
   end
   A(i, i) = 1 + A(i, i);
end
% розв'язок СЛАР    A*Y=F
Y = A\F;
return