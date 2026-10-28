function [y, alfa]  =  m_progn(a, c, b, f)
% метод немонотонної прогонки для розв'язання СЛАУ
%       з тридiагональною матрицею
% (метод Гауса з вибором провiдного елемента по рядку)
%    a(k)*y(k-1) - c(k)*y(k) + b(k)*y(k+1) = - f(k)
%         k= 1,2,...,n;  a(1)=b(n)=0 
%  вхiднi  данi:
%    a - пiддiагональ;
%   -c - дiагональ матрицi;
%    b - наддiагональ;
%   -f - права частина.
%  вихiднi данi:
%    y - розв'язок;
% alfa - прогоночнi коефiцiєнти.

n = length(c); n1 = n - 1;
% пряма прогонка:
alfa = zeros(1, n); y = alfa;
beta = alfa; ka = ones(1, n); te = alfa;
C = c(1); A = a(2); F = f(1); Q = f(2);
for k = 1 : n1
  j = k + 1;
  if abs(C) >= abs(b(k))
    alfa(j) = b(k) ./ C; beta(j) = F ./ C;
    C = c(j) - A .* alfa(j); F = Q + A .* beta(j);
    te(j) = ka(k); ka(j) = j;
    if k ~= n1
      j = j + 1; A = a(j); Q = f(j);
    end
  else
    alfa(j) = C ./ b(k); beta(j) = -F ./ b(k);
    C = c(j) .* alfa(j) - A; F = Q - c(j) .* beta(j);
    te(j) = j; ka(j) = ka(k);
    if k ~= n1
      j1 = j; j = j + 1;
      A = a(j) .* alfa(j1); Q = f(j) + a(j) .* beta(j1);
    end
  end  
end
% зворотня прогонка:
y(ka(n)) =F ./ C;
for k=n1 : -1 : 1
  j = k + 1; y(te(j)) = alfa(j) .* y(ka(j)) + beta(j);
end
return