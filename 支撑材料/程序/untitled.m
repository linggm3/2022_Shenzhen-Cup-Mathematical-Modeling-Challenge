clear;
clc;

figure;
syms x;
y = x^3 - 2*x + 2;
ax_1 = 0;
ax_2 = x * 9999999;
fplot(y);
hold on;
fplot(ax_1, 'k');
fplot(ax_2, 'k');
y_1 = 1 + (x-1);
y_2 = 2 - 2*x;
fplot(y_1, 'r');
fplot(y_2, 'r');
axis([-2, 2, -3, 5])
ylabel('Y')
xlabel('X')
title('y = x^3 - 2x + 2')

figure;
syms x;
y = tan(x);
x0 = pi / 2 - 0.02;
y_ = tan(x0) + (tan(x0)^2 + 1) * (x-x0);
fplot(y_, 'r');
ax_1 = 0;
ax_2 = x * 9999999;
hold on;
fplot(y, 'b');

fplot(ax_1, 'k');
fplot(ax_2, 'k');
axis([-4.7, 4.7, -10, 10])
ylabel('Y')
xlabel('X')
title('y = tan(x)')

figure;
syms x;
y = abs(x);
y_ = diff(y);
ax_1 = 0;
ax_2 = x * 9999999;
fplot(y);
hold on;
fplot(y_, 'r')
fplot(ax_1, 'k');
fplot(ax_2, 'k');
axis([-2, 2, -2, 3])
ylabel('Y')
xlabel('X')
title('y = |x|')

figure;
syms x;
y = sqrt(x);
x0 = 0.1
y_ = sqrt(x0) + 1/(2*x0^(1/2)) * (x-x0);
ax_1 = 0;
ax_2 = x * 9999999;
fplot(y);
hold on;
fplot(y_, 'r')
fplot(ax_1, 'k');
fplot(ax_2, 'k');
axis([-0.5, 1.5, -0.5, 1.5])
ylabel('Y')
xlabel('X')
title('y = sqrt(x)')

