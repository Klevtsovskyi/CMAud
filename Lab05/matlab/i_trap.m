function s = i_trap(a, b, f, ep, mki)
% «находженн€ з точнiстью ep значенн€ визначеного iнтеграла
%                      b
%                  s = S f(x) dx
%                      a
% з допомогою  ‘ методу трапецiй.
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
k = 0; h = B - A;
s0 = 0.5.*(f(A) + f(B)); s1 = 0; s = h.*s0;
F = 0;
while ~F & k < mki
   k = k + 1; ho = h; h = 0.5.*h;
   x = A + h;
   while x < B
      s1 = s1 + f(x); x = x + ho;
   end
   ss = h.*(s0 + s1);
   F = abs(ss - s)./3 < ep; s = ss;
end
s = z.*s;
return