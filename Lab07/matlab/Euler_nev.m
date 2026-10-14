function [Y, exitf] = Euler_nev(F, X, U0, epsx, epsf)
% Розв'язок U'(x) = F(x,U(x)), x є X,
%           U(X(1))=U0,
% неявним м.Ейлера на фіксованій сітці X.
    nX = length(X); N = nX-1;
    [n, m] = size(U0);
    if n >= m
        Y = [U0,zeros(n,N)];
    else
        Y = [U0',zeros(m,N)];
    end
    opt = optimset('Display','off','TolX',epsx,'TolFun',epsf);
    for k = 1 : N
       kp = k+1; x = X(kp); h = x - X(k);
       z = Y(:,k); Y0 = z;
       [YY, FX, exitf, output] = fsolve(@(T)(T-z-h.*F(x,T)), Y0, opt);
       if exitf ~= 1
           [k, x, YY, FX, exitf, output.iterations]
           return
       end
       Y(:,kp) = YY;
    end
end