function view_sol(Y, txt)
global X n
m = Y(1); M = m;
for i = 2:n
    de = Y(i);
    if de < m
        m = de;
    elseif de > M
        M = de;
    end
end
de = 1e-1; de = [X(1)-de, X(n)+de, m-de, M+de];
plot(X,Y,'b','LineWidth',2);
axis(de); grid on; title(txt);
xlabel('x'), ylabel('Y');
end