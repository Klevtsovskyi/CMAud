function Y  =  m_iim(Kfun, pfun, ffun, D)
% Наближений розв'язок крайової задачі
%   (K(x)*u'(x))' - p(x)*u(x) = -f(x),
%   K(x)>=K0>0, p(x)>=0, x є (a,b),
%  K(x), p(x), f(x) можуть бути розривні
% (-D(1)K(x)*u'(x)+D(2)u(x))|x=a = D(3),
%  (D(4)K(x)*u'(x)+D(5)u(x))|x=b = D(6).
% РС побудована інтегро-інтерполяційним методом.
% Вхiднi данi :
%   Kfun - функцiя з описом K(X(k));
%   pfun - функцiя з описом p(X(k));
%   ffun - функцiя з описом f(X(k));
%     D  - масив значень коефiцiєнтiв КУ:
%  D(1),D(4)=0 & D(2),D(5)=1  - КУ І роду,
%  D(1),D(4)=1 & D(2),D(5)>=0 - КУ IІ,III роду.
% Вихiднi данi :
%     Y  - масив значень наближеного розв'язку.
if nargin < 4
    error('M_IIM:NotEnoughInputs',...
          'Not enough input arguments.  See M_IIM.');
    return
end
if length(D) < 6
    error('M_IIM:NotEnoughInputs',...
          'Not enough input arguments D.  See M_IIM.');
    return
end
if D(1) ~= 0 & D(1) ~= 1
    error('M_IIM:NotEnoughInputs',...
          'Not enough input arguments D(1).  See M_IIM.');
    return
end
if D(4) ~= 0 & D(4) ~= 1
    error('M_IIM:NotEnoughInputs',...
          'Not enough input arguments D(4).  See M_IIM.');
    return
end
global X n h h05
hh = h.*h;  % допомiжнi змiннi
% формування коефiцiєнтiв
% A(k),C(k),B(k),F(k), k=1,...,n для СЛАР
A = zeros(1, n); B = A; C = A; F = A;
akp = Kfun(2);
if D(1)
    C(1) = akp + h.*(D(2) + h05.*pfun(1));
    B(1) = akp; F(1) = h.*(D(3) + h05.*ffun(1));
else
    C(1) = 1; F(1) = D(3);
end
for k = 2 : n-1
  ak = akp ; akp = Kfun(k+1);
  A(k) = ak; B(k) = akp;
  C(k) = ak + akp + hh.*pfun(k);
  F(k) = hh.*ffun(k);
end
if D(4)
    C(n) = akp + h.*(D(5) + h05.*pfun(n));
    A(n) = akp; F(n) = h.*(D(6) + h05.*ffun(n));
else
    C(n) = 1; F(n) = D(6);
end
% знах-ня Y з розв'язку СЛАР м.монотонної(правої) прогонки
[Y AL er] = m_progm(A, C, B, F);
return