function [X, Y] = m_grid2(pfun, ffun, a, b, D, n)
% Наближений розв'язок лiнiйного диференцiального рiвняння
%    u"(x) - p(x)*u(x) = - f(x),
%        p(x)>=0,  x є (a,b),
% для якого виконуються крайові умови 
%    D(1)*u'(a) + D(2)*u(a) = -D(3),
%    D(4)*u'(b) + D(5)*u(b) =  D(6).
%  (різницева схемах 2-го порядку точності).
% Для параметрiв D(1),D(2) i D(4),D(5) виконуються умова :
% abs(D(1)) + abs(D(2)) <> 0 & abs(D(4)) + abs(D(5)) <> 0.
% Вхiднi данi :
%   pfun   - функцiя з описом p(x);
%   ffun   - функцiя з описом f(x);
%   a, b   - границі вiдрiзку [a,b];
%     D    - масив коефiцiєнтiв крайових умов;
%     n    - число точок сiтки рiзницевої схеми.
% Вихiднi данi :
%     X    - масив значень вузлiв сiтки;
%     Y    - масив значень наближеного розв'язку.
if nargin<6
  error('M_GRID2:NotEnoughInputs',...
        'Not enough input arguments.  See M_GRID2.');
  return
end
if length(D)<6
  error('M_GRID2:NotEnoughInputs',...
        'Not enough input arguments D.  See M_GRID2.');
  return
end
if abs(D(1)) + abs(D(2)) == 0 | abs(D(4)) + abs(D(5)) == 0
  error('M_GRID2:NotEnoughInputs',...
        'Not enough input arguments D.  See M_GRID2.');
  return
end
if a >= b
  error('M_GRID2:error a>=b.','a>=b.  See M_GRID2.');
  return
end
% формування вузлiв
X = linspace(a, b, n); h = X(2) - X(1);
% допомiжнi змiннi
hh = h.*h; h2 = 2.*h;
% формування коефiцiєнтiв
%     A(k),C(k),B(k),F(k), k=1,...,n для СЛАР
A = ones(n, 1); B = A; C = A; F = A;
A(1) = 0; B(1) = 2.*D(1);  x = D(1).*hh;
C(1) = B(1) - D(2).*h2 + x.*pfun(a);
F(1) = D(3).*h2 + x.*ffun(a);
for k = 2 : n-1
  x = X(k); C(k) = 2 + hh.*pfun(x); F(k) = hh.*ffun(x);
end
A(n) = 2.*D(4); B(n) = 0; x = D(4).*hh;
C(n) = A(n) + D(5).*h2 + x.*pfun(b);
F(n) = D(6).*h2 + x.*ffun(b);
% знаходження сiткової функцiї при розв'язаннi СЛАР
% методом монотонної прогонки, а у випадку
% невиконання умов стiйкостi - методом немонотонної прогонки
[Y AL er] = m_progm(A, C, B, F);
if er
  [Y AL] = m_progn(A, C, B, F);
end
return