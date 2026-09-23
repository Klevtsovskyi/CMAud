function Le04_2
% Використання сплайнів для числового диференціювання
    clear all, close all, clc, warning off
    OM = true;   % true|false    Octave|MATLAB
    if OM
        pkg load symbolic   % для Octave
    end
    % означення неперервних u(x), u'(x), u''(x), u'''(x)
    u   = @(x)(sin(x));
    du  = @(x)(cos(x));
    d2u = @(x)(-sin(x));
    d3u = @(x)(-cos(x));

    % завдання сплайнової сітки X
    n = 16; X = linspace(0,6,n);

    % обчислення значень функції та її похідних на сітці X
    U = u(X); U1 = du(X); U2 = d2u(X); U3 = d3u(X);

    % обчислення кубічного сплайна (КС) з КУ IV-го роду
    % S31 - pp-форма КС в поліноміальній формі на кожному відрізку
    S31 = spline(X, U);

    % відображення S31 кусково-поліноміальної форми
    disp('pp-форма КС'); disp(S31);
    fprintf('%15s %36s %-s\n', 'Інтервал за х ',...
            'Коефіцієнти КС      ', '      Поліном');
    for i = 1:S31.pieces
        ci = sprintf('[%6.3f,%6.3f]', S31.breaks(i), S31.breaks(i+1));
        P = S31.coefs(i,:)';
        ck = sprintf(' %8.4f', P);
        fprintf('%s %s %s\n', ci, ck, poly_str(P));
    end
    clear ci P ck

    % сітка інтерполяції
    N = 10.*n+1; XX = linspace(0,6,N);

    % графіки: u(X),S31(XX;u); похибки R(XX)=S31(XX;u)-u(XX)
    frm = '%s %12.10f\n'; fprintf('\n');
    subplot(1,2,1);
    Y = ppval(S31,XX);
    plot(X,U,'k :',XX,Y,'b');
    grid on; xlabel('x'); ylabel('u');
    title('u(X) і S31(XX;u)');
    subplot(1,2,2);
    R = Y - u(XX); r = norm(R, inf);
    fprintf(frm, '||S31(x;u)  -u  ||c=', r);
    plot(XX,R);
    grid on; xlabel('x'); ylabel('R(x)');
    title('Похибка R(X) для u'); pause(7); close all

    % графіки: u'(X),S31'(XX;u); похибки R(XX)=S31'(XX;u)-u'(XX)
    subplot(1,2,1);
    if OM
        D = ppder(S31);
    else
        D = fnder(S31);   % 1-ша похідна від S31(x;u)
    end
    Y = ppval(D,XX);
    plot(X,U1,'k :',XX,Y,'b');
    grid on; xlabel('x'); ylabel('u''');
    title('u''(X) і S31''(XX;u)');
    subplot(1,2,2);
    R = Y - du(XX); r = norm(R, inf);
    fprintf(frm, '||S31''(x;u) -u'' ||c=', r);
    plot(XX,R);
    grid on; xlabel('x'); ylabel('R(x)');
    title('Похибка R(X) для u'''); pause(7); close all

    % графіки: u"(X),S31"(XX;u); похибки R(XX)=S31"(XX;u)-u"(XX)
    subplot(1,2,1);
    if OM
        D = ppder(D);
    else
        D = fnder(D);   % 2-га похідна від S31(x;u)
    end
    Y = ppval(D,XX);
    plot(X,U2,'k :',XX,Y,'b');
    grid on; xlabel('x'); ylabel('u''');
    title('u"(X) і S31"(XX;u)');
    subplot(1,2,2);
    R = Y - d2u(XX); r = norm(R, inf);
    fprintf(frm, '||S31"(x;u) -u" ||c=', r);
    plot(XX,R);
    grid on; xlabel('x'); ylabel('R(x)');
    title('Похибка R(X) для u"'); pause(7); close all

    % графіки: u"'(X),S31"'(XX;u); похибки R(XX)=S31"'(XX;u)-u"'(XX)
    subplot(1,2,1);
    if OM
        D = ppder(D);
    else
        D = fnder(D);   % 3-тя похідна від S31(x;u)
    end
    Y = ppval(D,XX);
    plot(X,U3,'k :',XX,Y,'b');
    grid on; xlabel('x'); ylabel('u"''');
    title('u"''(X) і S31"''(XX;u)');
    subplot(1,2,2);
    R = Y - d3u(XX); r = norm(R, inf);
    fprintf(frm, '||S31"''(x;u)-u"''||c=', r);
    plot(XX,R);
    grid on; xlabel('x'); ylabel('R(x)');
    title('Похибка R(X) для u"'''); pause(7); close all
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
