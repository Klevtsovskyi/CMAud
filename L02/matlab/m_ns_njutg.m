function [X, k, err] = m_ns_njutg(f_J, X, e, mki)
% Знаходження розв'язку нелiнiйної системи рiвнянь F(X) = 0
% iтерацiйним м.Ньютона + одновимірний пошук.
% Умова закiнчення iтерацiй : norm(F(X(k+1)))<e.
% вхiднi данi :
%   f_J    - функцiя, яка обчислює пару {F(X),J(X)} при fun=1,
%              а при fun=0 - {F(X)};
%   X      - початкове наближення до розв'язку;
%   e      - точнiсть (0<e<1);
%   mki    - максимальна кiлькiсть крокiв iтерацiй.
% вихiднi данi:
%   X      - фiнiшне наближення до розв'язку;
%   k      - кiлькiсть виконаних iтерацiй;
%   err    - 0=знайдено розв'язок, 1= k>mki.
global fun
if e < 0 | e > 1
    e = 1e-6;
end
for k = 1 : mki
   fun = 1; Xc = X; [F, J] = f_J(X); P = -J\F;
   fun = 0; v = norm(F);
   for i = 1: mki
       X = Xc + P; w = norm(f_J(X));
       if w <= v
           break
       end
       P = 0.5.*P;
   end
   if w < e
       err = 0; return
   end
end
err = 1;
return