function [x, k, err] = m_1simpit(ffun, x, e, mki)
% М.простої ітерації розв'язання рівняння  x = f(x)
% вхiднi данi :
%   ffun   - функцiя f(x);
%   x      - 0-наближення до розв'язку;
%   mki    - максимальна кiлькiсть крокiв iтерацiй;
%   e      - точнiсть (0<e<1);
% вихiднi данi:
%   x      - фiнiшне наближення до розв'язку;
%   k      - кiлькiсть виконаних iтерацiй.
%   err    - 0=знайдено, 1= k>mki.

if e < 0 | e > 1
    e = 1e-6;
end 
for k = 1 : mki
    r = ffun(x);
    if abs(r - x) < e
        err = 0; return
    end
    x = r;
end
err = 1;
return