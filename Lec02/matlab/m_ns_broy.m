function [X, k, err] = m_ns_broy(f_J, X, e, mki)
% Знаходження розв'язку нелiнiйної системи рiвнянь F(X) = 0
% iтерацiйним першим м.Бройдена, де X(0) - вiдоме.
% Умова закiнчення iтерацiй : norm(X(k+1)-X(k))<e.
% вхiднi данi :
%   f_J    - функцiя, яка при fun=0 обчислює {F(X)},
%            а при fun=1 - обчислює пару {F(X),J(X)},
%            де J(X) - якобіан (dFi(X)/dXj), i,j=1,...,n.
%   X      - початкове наближення до розв'язку;
%   e      - точнiсть (0<e<1);
%   mki    - максимальна кiлькiсть крокiв iтерацiй;
% вихiднi данi:
%   X      - фiнiшне наближення до розв'язку;
%   k      - кiлькiсть виконаних iтерацiй;
%   err    - 0=знайдено розв'язок, 1= k>mki.

global fun
fun = 1;
if e < 0 | e > 1
    e = 1e-6;
end
for k = 1 : mki
   if k == 1
       [O, J] = f_J(X); fun = 0;
   else
       O = f_J(X); J = J + O * Dx' ./dot(Dx,Dx);       
   end
   Dx = -J\O; X = X + Dx;
   if norm(Dx) < e
       err = 0; return
   end
end
err = 1;
return