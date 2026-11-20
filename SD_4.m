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

F = @(x, y) (1/epsilon) * ( -x + y ./ (1 + y) ) - y ./ ((1 + y).^2);
G = @(x, y) -y; 

[X, Y] = meshgrid(linspace(0, 1.5, 16));

U = F(X, Y);
V = G(X, Y);

y_min = min(ya15s(:)) - 0.05;
y_max = max(ya15s(:)) + 0.05;

figure; 
quiver(X, Y, U, V)
grid on
hold on
plot(ya15s(:, 1), ya15s(:, 2), LineWidth=2)
fimplicit(F, [0 1.5 0 1.5], 'k', LineWidth=1)
xline(0, 'k')
axis([0 1.5 0 1.5])
title('Spazio delle fasi', Interpreter='latex')
xlabel('y1')
ylabel('y2')


figure; 
tiledlayout(1, 2)

nexttile
h1 = animatedline('LineWidth', 1.5, Color='blue');
hold on 
h2 = animatedline('LineWidth', 1.5, Color='red');
axis([0, T, y_min y_max])
grid on
title(sprintf('Evoluzione temporale, $\\epsilon$ = %g', epsilon), Interpreter='latex')

nexttile
h3 = animatedline('LineWidth', 1.5, Color='red');
hold on
quiver(X, Y, U, V)
fimplicit(F, [0 1.5 0 1.5], 'k', LineWidth=1)
yline(0, 'k', LineWidth=1)
axis([ y_min y_max y_min y_max])
grid on
title(sprintf('Spazio delle fasi, $\\epsilon$ = %g', epsilon), Interpreter='latex')
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


%%
clear; clc;

epsilon_test = 10.^(-1:-1:-5);

y01 = 1.5; 
y02 = 1; 
T = 10; 

ITER = zeros(length(epsilon_test), 3);
ERR = zeros(length(epsilon_test), 3);
TIME = zeros(length(epsilon_test), 3); 

for idx = 1:length(epsilon_test)
    epsilon = epsilon_test(idx);

    f = @(t, y) [(1/epsilon) * ( -y(1) + y(2) ./ (1+y(2)) ) - y(2) ./ ((1+y(2)).^2); ...
        -y(2)];

    y1 = @(t)  ( y01 - y02/(1+y02) ) * exp(-t/epsilon) + y02*exp(-t)./(1+y02*exp(-t));
    y2 = @(t) y02*exp(-t);

    y = @(t) [y1(t) y2(t)];

    tspan = [0, T];
    y0 = [y01 y02];

    tStart113 = cputime;
    [t113, ya113] = ode113(f, tspan, y0);
    TIME(idx, 1) = max(cputime - tStart113, 1e-2); % uso il max solo per ragioni di plot 

    tStart15s = cputime;
    [t15s, ya15s] = ode15s(f, tspan, y0);
    TIME(idx, 2) = max(cputime - tStart15s, 1e-2);

    tStart45 = cputime;
    [t45, ya45] = ode45(f, tspan, y0);
    TIME(idx, 3) = max(cputime - tStart45, 1e-2);

    err113 = norm(sqrt(sum((y(t113)-ya113).^2, 2))) / norm(sqrt(sum(y(t113).^2, 2)));
    err15s = norm(sqrt(sum((y(t15s)-ya15s).^2, 2))) / norm(sqrt(sum(y(t15s).^2, 2)));
    err45 = norm(sqrt(sum((y(t45)-ya45).^2, 2))) / norm(sqrt(sum(y(t45).^2, 2)));

    ITER(idx, 1) = length(t113);
    ITER(idx, 2) = length(t15s);
    ITER(idx, 3) = length(t45);

    ERR(idx, 1) = err113;
    ERR(idx, 2) = err15s;
    ERR(idx, 3) = err45;

end

figure; 
loglog(epsilon_test, ERR, '--o', LineWidth=2)
grid on
legend('ode113', 'ode15s', 'ode45', Location='best')
xlabel('$\epsilon$', Interpreter='latex')
title('Errore relativo', Interpreter='latex')

figure; 
loglog(epsilon_test, ITER, '--o', LineWidth=2)
grid on
legend('ode113', 'ode15s', 'ode45', Location='best')
xlabel('$\epsilon$', Interpreter='latex')
title('Numero passi temporali', Interpreter='latex')

figure; 
loglog(epsilon_test, TIME, '--o', LineWidth=2)
grid on
legend('ode113', 'ode15s', 'ode45', Location='best')
xlabel('$\epsilon$', Interpreter='latex')
title('Tempo di cpu', Interpreter='latex')










