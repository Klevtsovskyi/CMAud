function view_solu(Y, U, ni, txt, s)
global X n h
m = U(1); M = m;
for i = 2:n
    de = U(i);
    if de < m
        m = de;
    elseif de > M
        M = de;
    end
end
i = fix(linspace(1, n, ni));
de = 1e-1; de = [X(1)-de, X(n)+de, m-de, M+de];
plot(X,Y,'b',X(i),U(i),'k*','LineWidth',2);
axis(de); grid on; title(txt);
xlabel('x'), ylabel('Y,u');
legend('Y','u','Location','Best');
if s, pause(s)
else, pause
end
close all
U = abs(Y-U);
fprintf('h^2=%11.3e Похибка в нормі C=%13.5e\n',h.*h,max(U));
plot(X,U,'r','LineWidth',1.5); grid on;
title('Абсолютна похибка');
xlabel('x'), ylabel('|Y-u(X)|');
end