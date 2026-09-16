function [x, k, err] = m_1dihot(ffun, D, e, mki)
% Знаходження кореня нелiнiйного рiвняння  f(x) = 0
% методом діхотомії.
% вхiднi данi :
%   D      - діапазон кореня;
%   e      - точнiсть (0<e<1);
%   mki    - максимальна кiлькiсть крокiв iтерацiй;
%   ffun   - функцiя f(x).
% вихiднi данi:
%   x      - фiнiшне наближення до кореня;
%   k      - кiлькiсть виконаних iтерацiй;
%   err    - -1= не виконується умова D(1)<D(2)&f(D(1))*f(D(2))<=0;
%             0= знайдено,
%             1= k>mki.

if e < 0 | e > 1
    e = 1e-6;
end
l = D(1); r = D(2); fl = ffun(l); fr = ffun(r);
if fl.*fr <= 0
    err = 1;
    for k = 1 : mki
        x = 0.5.*(l+r); fx = ffun(x);
        if fl.*fx <= 0
            r = x; fr = fx;
        else
            l = x; fl = fx;
        end
        if abs(r-l) < e
            err = 0; break
        end
    end
    x = 0.5.*(l+r);
else
    err = -1; k = 0;
end
return