function Pr04_2
   clear all, close all, clc, warning off
   
   a = -1; b = 1;
   N = [5, 11, 17];
   C = {'b -','g -','r -'};
   
   T = 'Пр. С.Н.Бернштейна';
   u = @(t)(abs(t));
   viewIP(u, a, b, N, C, T);
   pause; close all; clear u

   T = 'Приклад К.Рунге';
   u = @(t)(1./(1+25.*t.^2));
   viewIP(u, a, b, N, C, T);
end

function viewIP(u, a, b, N, C, T)
   nN = length(N); s = 5; 
   for k = 1 : nN
      n = N(k); nm = n-1; Ns = nm*s+1; 
      X = linspace(a, b, n);    % сітка інтерполяційна
      x = linspace(a, b, Ns);   % сітка для обчислень ("згущена")
      U = u(X);                 % значення u(t) на сітці X
      P = polyfit(X, U, nm);    % побудова інтерп.полінома P на X
      Y = polyval(P, x);        % значення ІП P в точках x

      subplot(1,3,1);
      plot(x,Y,C{k}); hold on;
      if k == nN
         plot(X,U,'k. '); hold on;
         grid on; axis equal;
         xlabel('X,x'); ylabel('U,Y'); title(T);
      end

      subplot(1,3,2);
      S11 = interp1(X, U, x);
      plot(x,S11,C{k}); hold on;
      if k == nN
         plot(X,U,'k. '); hold on;
         grid on; axis equal;
         xlabel('X,x'); ylabel('U,S11'); title(T);
      end

      subplot(1,3,3);
      S31 = spline(X, U, x);
      plot(x,S31,C{k}); hold on;
      if k == nN
         plot(X,U,'k. '); hold on;
         grid on; axis equal;
         xlabel('X,x'); ylabel('U,S31'); title(T);
      end
   end
   hold off;
end