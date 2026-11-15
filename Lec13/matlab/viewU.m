function viewU(s, X, Y, U, kc, UU)
   [Mx, My] = meshgrid(X, Y); close all;
   mesh(Mx, My, U);
   e2 = 'Різницевий розв’язок';
   title(e2); xlabel('x'); ylabel('y'); zlabel('U');
   shading interp; colormap(pink); %view(40,20);
   if s, pause(s);
   else, pause;
   end
   close all;
   if kc
      contour(Mx, My, U, kc);
      xlabel('x'); ylabel('y'); title(e2); colormap hsv;
      if s, pause(s);
      else, pause;
      end
      close all;
   end

   if nargin == 5
       return
   end

   mesh(Mx, My, UU);
   e2 = 'Точний розв’язок';
   title(e2); xlabel('x'); ylabel('y'); zlabel('u');
   shading interp; colormap(pink); %view(40,20);
   if s, pause(s);
   else, pause;
   end
   close all;
   if kc
      contour(Mx, My, UU, kc);
      xlabel('x'); ylabel('y'); title(e2); colormap hsv;
      if s, pause(s);
      else, pause;
      end
      close all;
   end

   Uu = U - UU;
   mesh(Mx, My, Uu);
   e2 = 'Похибка U-u';
   title(e2); xlabel('x'); ylabel('y'); zlabel('U-u');
   shading interp; colormap(pink);  %view(40,20);
   if s, pause(s);
   else, pause;
   end
end