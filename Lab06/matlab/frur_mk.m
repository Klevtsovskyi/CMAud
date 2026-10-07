function [X Y] = frur_mk(a, b, nkf, n, k_dkfun, f_dffun, Y, e, mki, m)
% Розв'язок нелiнiйного iнтегрального рiвняння
% типу Фредгольма-Урисона методом замiни
% iнтеграла скiнченною сумою
%     b
%     S K(x,s,u(s))ds = f(x,u(x)),   x є [a,b].
%     a
% вхiднi данi :
%     a,b - нижня і верхня границя iнтеграла;
%     nkf - 0 = КФ методу трапецiй,
%           1 = КФ методу Сiмпсона;
%       n - число точок сiтки. Повинно задовольняти умовi:
%            n > 1,                 при nkf=0,
%            n > 2 & mod(n,2)==1,   при nkf=1;
% k_dkfun - функцiя з описом ядра K(x,s,y) i dK(x,s,u)/du;
% f_dffun - функцiя з описом правої частини f(x,y) i df(x,u)/du;
%       Y - початкове наближення до розв'язку;
%       e - точнiсть (0<e<1);
%     mki - максимальна кiлькiсть крокiв iтерацiй;
%       m - метод розв'язку нелінійної системи:
%            m = 1 - градієнта;
%            m = 2 - Ньютона;
%            m = 3 - перші 4 ітерації по градієнтному, далі по Ньютона;
% вихiднi данi:
%  X    - масив вузлiв сiтки;
%  Y    - масив наближеного розв'язку.
if nargin<10
   error('FRUR_MK:NotEnoughInputs',...
         'Not enough input arguments.  See FRUR_MK.');
   return
end
if a >= b
   error('FRUR_MK:error a >= b.',...
          'a >= b.  See FRUR_MK.');
   return
end
switch nkf
   case 0
      if n < 2
         error('FRUR_MK:error nkf = 0 & n < 2.',...
                 'nkf = 0 & n < 2.  See FRUR_MK.');
         return
      end
   case 1
      if n < 3
         error('FRUR_MK:error nkf = 1 & n < 3.',...
                 'nkf = 1 & n < 3.  See FRUR_MK.');
         return
      else
         if mod(n,2) == 0
            error('FRUR_MK:error nkf = 1 & mod(n,2) == 0.',...
                    'nkf = 1 & mod(n,2) == 0.  See FRUR_MK.');
            return
         end
      end
   otherwise
      error('FRUR_MK:Unknown method.',...
            'Unknown method.  See FRUR_MK.');
      return
end
if m ~= 1 & m ~= 2 & m ~= 3
   error('FRUR_MK:error m ~=[1,2,3].',...
          'Unknown method m (m ~є [1,2,3]).  See FRUR_MK.');
   return
end
% формування вузлiв сітки
X = linspace(a,b,n); h = X(2)-X(1); n1 = n-1;
d = Y; F = Y; J = zeros(n,n);
% формування коефiцiєнтiв P(i), i=1,...,n КФ
P = Y;
switch nkf
   case 0
      P(1) = 0.5.*h;
      for i=2:n1
         P(i) = h;
      end
      P(n) = P(1);
   case 1
      P(1) = h./3; c2 = 2.*P(1);  c1 = 2.*c2;
      for i=2:2:n1
         P(i) = c1; P(i+1) = c2;
      end
      P(n) = P(1);
end
clear c1 c2 h
s = 0; Fl = 1;
while Fl & s < mki
% знаходження розв'язку нелінійної системи :
% формування матрицi Якобі J i вектора F
   for i = 1:n
      xi = X(i); sk = 0;
      for j = 1:n
         p = P(j);
         [k,dk] = k_dkfun(xi,X(j),Y(j));
         J(i,j) = p.*dk; sk = sk + p.*k;
      end
      [f,df] = f_dffun(xi,Y(i));
      J(i,i) = J(i,i) - df;
      F(i) = sk - f;
   end
   if m == 1 | m == 3 & s < 3
% метод градієнта
      d = J' * F; v = J * d;
      sk = dot(F,v)./dot(v,v);
      d = sk.*d; 
   else       
% метод Ньютона
      d = J\F;
   end
   Y = Y - d; s = s + 1
   sk = norm(d); Fl = sk > e;
end
clear d J F sk P xi p k dk f df
Y = Y';
return