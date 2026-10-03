function Pr06_1
   clear all, close all, clc

   n = 7; a = 0; b = 1; ng = 4;
   L = b-a; h = L./(n-1);
   I = int(sym('sin(pi*x)'),sym('x'),a,b); I = eval(I);
   % Обчислення КФ :
   Cr = i_crectn(a, b, @f, n); % центр.прямокутників
   T  = i_trapn(a, b, @f, n);  % трапецій
   S  = i_simpn(a, b, @f, n);  % Сімпсона
   G  = i_gauss(a, b, @f, ng); % Гаусса

   h2 = h.*h; h4 = h2.*h2;
   M2 = pi.*pi; M4 = M2.*M2; p = L.*h2.*M2./12;
   % Друк таблиці результатів
   fprintf('%12s%19.12f\n','Точне знач.=',I);
   fprintf('%12s %12s %12s %12s\n',...
           'Назва КФ','Наближене значення',...
           'Абсолютна  похибка','Теор.оцінк.похибки');
   fr = '%12s%19.12f%19.10e%19.10e\n';
   NP = sprintf('n=%u',n); NG = sprintf('n=%u',ng);
   fprintf(fr,['центр.пр ',NP],Cr,abs(I-Cr),p./2);        
   fprintf(fr,['трапецій ',NP], T,abs(I-T),p);
   fprintf(fr,['Сімпсона ',NP], S,abs(I-S),L.*h4.*M4./180);
   fprintf(fr,['Гаусса   ',NG], G,abs(I-G),1e-8.*M4.*M4);
end

function y = f(x)
    y = sin(pi.*x);
end