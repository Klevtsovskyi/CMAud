function s = i_simp(a, b, f, ep, mki)
% Знаходження з точнiстью ep значення визначеного iнтеграла
%                      b
%                  s = S f(x) dx
%                      a
% з допомогою КФ методу Сiмпсона.
% Вхiднi данi :
%    a - нижня границя iнтеграла;
%    b - верхня границя iнтеграла;
%    f - функцiя з описом f(x);
%  mki - максимальна кiлькiсть згущення сiтки.
% Вихiднi данi :
%    s - значення iнтегралу.
A = a; B = b; z = 1;
if a > b
   A = b; B = a; z = -1;
end
k = 0; h = B - A; h = 0.5.*h;
s1 = f(A) + f(B); s2 = 0; s4 = f(0.5.*(B+A));
s = h.*(s1 + 4.*s4)./3;
F = 0;
while ~F & k < mki
    k = k+1; ho = h; h = 0.5.*h;
    s2 = s2 + s4;
    s4 = 0; x = A + h;
    while x < B
        s4 = s4 + f(x); x = x + ho;
    end
    ss = h.*(s1 + 4.*s4 + 2.*s2)./3;
    F = abs(ss - s)./15 < ep; s = ss;
end
s = z.*s;
return