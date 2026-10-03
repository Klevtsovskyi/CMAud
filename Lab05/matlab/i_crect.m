function s = i_crect(a, b, f, ep, mki)
% «находженн€ з точнiстью ep значенн€ визначеного iнтеграла
%                      b
%                  s = S f(x) dx
%                      a
% з допомогою  ‘ методу центр.пр€мокутник≥в.
% ¬хiднi данi :
%  a, b Ц пром≥жок iнтегруванн€;
%     f Ц функцi€ з описом f(x);
%    ep Ц точнiсть (0 < ep < 1);
%   mki Ц максимальна кiлькiсть згущенн€ сiтки.
% ¬ихiднi данi :
%     s Ц значенн€ iнтегралу.
A = a; B = b; z = 1;
if a > b
   A = b; B = a; z = -1;
end
k = 0; h = B - A; s = h.*f(0.5.*(A+B));
F = 0;
while ~F & k < mki
   k = k + 1; h = 0.5.*h;
   x = A + 0.5.*h; s1 = 0;
   while x < B
      s1 = s1 + f(x); x = x + h;
   end
   ss = h.*s1;
   F = abs(ss - s)./3 < ep; s = ss;
end
s = z.*s;
return