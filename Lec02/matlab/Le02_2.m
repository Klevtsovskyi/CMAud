%%Le02_2.m
% Знаходження кореня нелінійного рівняння f(x)=0.
close all, clear all, clc, warning off

f = @(x)(x.*x-2.*x.*sin(x)-1);
F = @(x)(f(x).^2);

ezplot(F,[-3,3]); hold on; ezplot(f,[-3,3]); grid on

D = input('Інтервал, де знаходиться корінь=');
z = input('Наближення=');

x = fminbnd(F,D(1),D(2));
fprintf('Мінімізація  : F(%.6f)=%.6e\n', x, F(x));

x = fzero(f, z);  % 1-й варіант виклику
fprintf('Функція fzero: f(%.6f)=%.6e\n', x, f(x));
%з розширеним переліком результатів і виведенням процесу пошуку
[x, fx, ef] = fzero(f, D, optimset('Display','iter'));
fprintf('f(%.6f)=%.6e Код завершення=%2d\n', x, fx, ef);