function Y = Euler_ev(F, X, U0)
% Розв'язок U'(x) = F(x,U(x)), x є X,
%           U(X(1))=U0,
% явним м.Ейлера на фіксованій сітці X.
    nX = length(X); N = nX-1;
    [n, m] = size(U0);
    if n >= m
        Y = [U0,zeros(n,N)];
    else
        Y = [U0',zeros(m,N)];
    end
    for k = 1 : N
       kp = k+1; x = X(k); h = X(kp) - x;
       z = Y(:,k); Y(:,kp) = z + h.*F(x,z);
    end
end