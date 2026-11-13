clear; clc; close all;

epsilon = 1e-1; 

y01 = 1.5; 
y02 = 1; 
T = 10; 

f = @(t, y) [(1/epsilon) * ( -y(1) + y(2) ./ (1+y(2)) ) - y(2) ./ ((1+y(2)).^2); ...
    -y(2)];

y1 = @(t)  ( y01 - y02/(1+y02) ) * exp(-t/epsilon) + y02*exp(-t)./(1+y02*exp(-t));
y2 = @(t) y02*exp(-t);

y = @(t) [y1(t) y2(t)];

tspan = [0, T];
y0 = [y01 y02];

[t113, ya113] = ode113(f, tspan, y0);
[t15s, ya15s] = ode15s(f, tspan, y0);
[t45, ya45] = ode45(f, tspan, y0);

err113 = norm(sqrt(sum((y(t113)-ya113).^2, 2))) / norm(sqrt(sum(y(t113).^2, 2)));
err15s = norm(sqrt(sum((y(t15s)-ya15s).^2, 2))) / norm(sqrt(sum(y(t15s).^2, 2)));
err45 = norm(sqrt(sum((y(t45)-ya45).^2, 2))) / norm(sqrt(sum(y(t45).^2, 2)));


% plot campo di direzioni 

F= @(x, y) (1/epsilon) * ( -x + y ./ (1 + y) ) - y ./ ((1 + y).^2);
G = @(x, y) -y; 

[X, Y] = meshgrid(linspace(0, 1.5, 40));

U = F(X, Y);
V = G(X, Y);

y_min = min(ya15s(:)) - 0.05;
y_max = max(ya15s(:)) + 0.05;

figure; 
quiver(X, Y, U, V)
grid on
hold on
plot(ya15s(:, 1), ya15s(:, 2))
axis([ y_min y_max y_min y_max])


figure; 
tiledlayout(1, 2)

nexttile
h1 = animatedline('LineWidth', 1.5, Color='blue');
h2 = animatedline('LineWidth', 1.5, Color='red');
axis([0, T, y_min y_max])
grid on
title('Evoluzione temporale', Interpreter='latex')

nexttile
h3 = animatedline('LineWidth', 1.5, Color='black');
axis([ y_min y_max y_min y_max])
grid on
title('Spazio delle fasi', Interpreter='latex')
xlabel('y1')
ylabel('y2')

for k = 1:length(t15s)
    addpoints(h1, t15s(k), y1(t15s(k)))
    addpoints(h2, t15s(k), y2(t15s(k)))

    addpoints(h3, y1(t15s(k)), y2(t15s(k)))

    drawnow
    pause(0.1)
end
drawnow
