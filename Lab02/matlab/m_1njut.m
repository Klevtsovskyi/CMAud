function [x, k, err] = m_1njut(ffun, x, e, mki, dffun, N)
% Знаходження кореня нелiнiйного рiвняння  f(x) = 0
% iтерацiйним методом Ньютона  x(n+1)=x(n)-f(x(n))/f'(x(n)),
% де значення x(0) вiдоме.
% Умова закiнчення iтерацiй : abs(x(n+1)-x(n))<e.
% вхiднi данi :
%   x      - початкове наближення до кореня;
%   e      - точнiсть (0<e<1);
%   mki    - максимальна кiлькiсть крокiв iтерацiй;
%   ffun   - функцiя f(x);
%   dffun  - функцiя f'(x);
%   N       - кратність кореня.
% вихiднi данi:
%   x      - фiнiшне наближення до кореня;
%   k      - кiлькiсть виконаних iтерацiй;
%   err    - 0=знайдено, 1= k>mki, 2= інші помилки.

if e < 0 | e > 1
    e = 1e-6;
end
if N < 1
    k = 0; err = 2; return
end
for k = 1 : mki
    o = N.*ffun(x)./dffun(x); x = x-o;
    if abs(o) < e
        err = 0; return
    end
end
err = 1; return