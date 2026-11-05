function [y, alfa, err]  =  m_progm(a, c, b, f)
% метод монотонної (правої) прогонки для розв'язання
%       СЛАР з тридiагональною матрицею
%    a(k)*y(k-1) - c(k)*y(k) + b(k)*y(k+1) = - f(k),
%                a(1) = b(n) = 0,
% при виконаннi умов :
%            |c(k)| >= |a(k)| + |b(k)|,
%                k= 1,2,...,n.
%  вхiднi  данi:
% a, -c, b - пiддiагональ, дiагональ і наддiагональ матрицi СЛАР;
%       -f - права частина.
%  вихiднi данi:
%        y - розв'язок;
%     alfa - прогоночнi коефiцiєнти;
%      err - 1= умови застосування методу не виконуються i
%             розв'язок не знайдено,
%            0= розв'язок знайдено.
n = length(c); err = 0; k = 0;
while k < n & ~err
   k = k + 1;
   err = err | abs(c(k)) < abs(a(k)) + abs(b(k));
end
alfa = zeros(1, n); y = alfa;
if err
   return
end
% пряма прогонка :
nm = n - 1; beta = alfa;
for k = 1 : nm
   j = k + 1; z = c(k) - a(k).*alfa(k);
   alfa(j) = b(k)./z; beta(j) = (f(k) + a(k).*beta(k))./z;
end
% зворотна прогонка :
y(n) = (f(n) + a(n).*beta(n))./(c(n) - a(n).*alfa(n));
for k = nm : -1 : 1
   j = k + 1; y(k) = alfa(j).*y(j) + beta(j);
end
return