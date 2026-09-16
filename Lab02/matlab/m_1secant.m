function [x, k, err] = m_1secant(ffun, x, e, mki, x0)
%М.січних знаходження розв'язку рiвняння f(x)=0
% вхiднi данi :
%   f      - функцiя f(x);
%   x      - 1-наближення до розв'язку;
%   x0     - 0-наближення до розв'язку;
%   mki    - максимальна кiлькiсть крокiв iтерацiй;
%   e      - точнiсть (0<e<1);
% вихiднi данi:
%   x      - фiнiшне наближення до розв'язку;
%   k      - кiлькiсть виконаних iтерацiй.
%   err    - 0=знайдено, 1= k>mki.

if e < 0 | e > 1
    e = 1e-6;
end
fxm = ffun(x0);
for k = 1 : mki
    fx = ffun(x);
    o = (x-x0)./(fx-fxm).*fx;
    x0 = x; x = x - o;
    if abs(o) <= e
        err = 0; return
    end
    fxm = fx;
end
err = 1;
return