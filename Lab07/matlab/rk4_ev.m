function Y = rk4_ev(F, X, U0)
% Розв'язок U'(x) = F(x,U(x)), x є X,
%           U(X(1))=U0,
% явним 4-стадійним м.Рунге-Кутта 4-го порядку
% на фіксованій сітці X.
    nX = length(X); N = nX-1;
    [n, m] = size(U0);
    if n > m
        Y = [U0,zeros(n,N)];
    else
        Y = [U0',zeros(m,N)];
    end
    a1 = 1./6; a2 = 1./3;
    for k = 1 : N
       kp = k+1; x = X(k); xp = X(kp); h = xp - x;
       g = 0.5.*h; z = Y(:,k); K1 = F(x,z);
       xt = x+g; K2 = F(xt,z+g.*K1);
       K3 = F(xt,z+g.*K2); K4 = F(xp,z+h.*K3);
       Y(:,kp) = z + h.*(a1.*K1+a2.*K2+a2.*K3+a1.*K4);
    end
end