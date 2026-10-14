function Y = rk2_ev(F, X, U0, gam)
% Розв'язок U'(x) = F(x,U(x)), x є X,
%           U(X(1))=U0,
% явним 2-стадійним м.Рунге-Кутта 2-го порядку
% на фіксованій сітці X.
% Сімейство схем задається 0 < gam <= 1.
    if gam <= 0 | gam > 1
        gam = 1;  % схеми предиктор-коректор
    end
    nX = length(X); N = nX-1;
    [n, m] = size(U0);
    if n > m
        Y = [U0,zeros(n,N)];
    else
        Y = [U0',zeros(m,N)];
    end
    A1 = 1-gam; A2 = gam; g2 = 0.5./gam;
    for k = 1 : N
       kp = k+1; x = X(k); h = X(kp) - x;
       hg = h.*g2; z = Y(:,k);
       K1 = F(x,z); K2 = F(x+hg,z+hg.*K1);
       Y(:,kp) = z + h.*(A1.*K1+A2.*K2);
    end
end