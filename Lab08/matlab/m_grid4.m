function [X, Y] = m_grid4(p, ffun, dffun, a, b, D, n)
% Наближений розв'язок лiнiйного диференцiального рiвняння
%    u"(x) - p*u(x) = - f(x),
%  p є Const, p>0,  x є (a,b),
% для якого виконуються крайові умови 
%    D(1)*u(a) + D(2)*u'(a) = -D(3),
%    D(4)*u(b) + D(5)*u'(b) =  D(6).
%  (різницева схемах 4-го порядку точності).
% Для параметрiв D(1),D(2) i D(4),D(5) виконуються умова :
% abs(D(1)) + abs(D(2)) <> 0 & abs(D(4)) + abs(D(4)) <> 0.
% Вхiднi данi :
%   p    - коєффiцієнт в рівнянні;
%  ffun  - функцiя обчислення [f(x),f"(x)];
% dffun  - функцiя обчислення f'(x);
%   a, b - границі вiдрiзку [a,b];
%   D    - масив значень коефiцiєнтiв крайових умов;
%   n    - число точок сiтки рiзницевої схеми.
% Вихiднi данi :
%   X      - масив значень вузлiв сiтки;
%   Y      - масив значень наближеного розв'язку.
if nargin<7
  error('M_GRID4:NotEnoughInputs',...
        'Not enough input arguments.  See M_GRID4.');
  return
end
if length(D)<6
  error('M_GRID4:NotEnoughInputs',...
        'Not enough input arguments D.  See M_GRID4.');
  return
end
if abs(D(1)) + abs(D(2)) == 0 | abs(D(4)) + abs(D(5)) == 0
  error('M_GRID4:NotEnoughInputs',...
        'Not enough input arguments D.  See M_GRID4.');
  return
end
if a >= b
  error('M_GRID4:error a>=b.','a>=b.  See M_GRID4.');
  return
end
% формування вузлiв
X = linspace(a, b, n); h = X(2) - X(1);
% допомiжнi змiннi
hh = h.*h; hh6 = hh./6; hh12 = hh./12;
h4 = 0.25.*h; h2 = 0.5.*h;
ph6 = 1 + p.*hh6; ph12 = 1 + p.*hh12;
e1 = ph12.*hh; hv = hh12.*hh; x = 2 + p.*e1;
% формування коефiцiєнтiв
%   A(k),C(k),B(k),F(k), k=1,...,n для СЛАР
A = ones(1, n); B = A; C = x.*A; F = A;
% враховуємо крайові умови при x = a
A(1) = 0; B(1) = D(1)./h;
C(1) = B(1) + h2.*p.*D(1).*ph12 - D(2).*ph6;
[f, f2] = ffun(a); f1 = dffun(a);
F(1) = D(3).*ph6 + D(1).*(h2.*ph12.*f + hh6.*(f1 + h4.*f2));
% обчислюємо праву частину для i=2:n-1 :
for i = 2 : n-1
   x = X(i); [f, f2] = ffun(x); F(i) = e1.*f + hv.*f2;
end
% враховуємо крайові умови при x = b
A(n) = D(4)./h; B(n) = 0;
C(n) = A(n) + h2.*p.*ph12.*D(4) + ph6.*D(5); 
[f, f2] = ffun(b); f1 = dffun(b);
F(n) = ph6.*D(6) + D(4).*(h2.*ph12.*f - hh6.*(f1 - h4.*f2));
% розв'язок СЛАР методом немонотонної прогонки :
[Y, alfa] = m_progn(A, C, B, F);
return