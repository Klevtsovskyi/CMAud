function Le05_1
% Використання сплайнів для числового інтегрування
    clear all, close all, clc, warning off
    OM = true; %  true|false Octave|MATLAB
    if OM
        pkg load symbolic   % для Octave
    end
    % означення u(x)
    u   = @(x)(sin(x));
    % завдання сплайнової сітки X
    n = 16; X = linspace(0,6,n);
    % обчислення значень функції
    U = u(X);
    % обчислення кубічного сплайна (КС) з КУ IV-го роду
    % S31 - pp-форма КС в поліноміальній формі на кожному відрізку
    S31 = spline(X, U);
    % обчислення невизначеного інтеграла
    if OM
        In = ppint(S31);
    else
        In = fnint(S31);
    end
    fprintf('%15s %36s %-s\n', 'Інтервал за х ',...
            'Коефіцієнти КС      ', '      Поліном');
    for i = 1:In.pieces
        ci = sprintf('[%6.3f,%6.3f]', In.breaks(i), In.breaks(i+1));
        P = In.coefs(i,:)';
        ck = sprintf(' %8.4f', P);
        fprintf('%s %s %s\n', ci, ck, poly_str(P));
    end
    clear ci P ck
    % символьне обчислення інтеграла
    syms x;
    f = sin(x);
    IR = int(f);
    Ir = subs(IR,x,6)-subs(IR,x,0);
    if OM
        Ir = eval(Ir);   % для  Octave  !
    end
    fprintf('\n%-15s %18.10e\n', 'Точне значення=', Ir);
    fprintf('%-15s %18s %-s\n', 'Функція',...
            'Наближене значення', 'Похибка');
    frm = '%-15s %18.10e %18.10e\n';
    if OM
        I = ppval(In,6)-ppval(In,0);  % використання сплайна
        %I = eval(I);   % для  Octave  !
    else
        I = fnval(In,6)-fnval(In,0);  % використання сплайна
    end
    fprintf(frm,'fnint(S31)', I, I-Ir);
    S = trapz(X, U);              % м.трапецій
    fprintf(frm,'trapz(X,U)', S, S-Ir);
    S = quad(u, 0, 6);            % м.Сімпсона
    fprintf(frm,'quad(u,0,6)', S, S-Ir);
    S = quadl(u, 0, 6);           % м.Сімпсона підвищеної точності
    fprintf(frm,'quadl(u,0,6)', S, S-Ir);
    % сітка інтерполяції
    N = 10.*n+1; XX = linspace(0,6,N);
    % графічне порівняння
    subplot(1,2,1);
    Y = ppval(In,XX);
    for i=1:N
        C = subs(IR,XX(i));
        if OM
            C = eval(C);   % для  Octave  !
        end
        U(i) = C;
    end
    C = Y(1)-U(1); Y = Y-C;
    plot(XX,U,'k :',XX,Y,'b');
    grid on; xlabel('x'); ylabel('int(u)');
    title('Невизначений інтеграл');
    subplot(1,2,2);
    R = Y-U; r = norm(R, inf);
    fprintf('\n%s %12.10f\n', '||In-int(u)||c=', r);
    plot(XX,R);
    grid on; xlabel('x'); ylabel('R(x)');
    title('Похибка R(X) для int(u)');
end

function S = poly_str(P)
% Побудова рядка S по поліному P
    syms x;
    n = length(P); k = n; S = 0;
    for i=1:n
        k = k-1; S = S + P(i)*x^k;
    end
    S = vpa(simplify(S),4); S = char(S);
end
