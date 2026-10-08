clear; clc; close all; 

%% punto a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulare il sistema dinamico dell'esercizio e plottare le soluzioni y1,
% y2, r=sqrt(y1^2 + y2^2) e nello spazio delle fasi (y1, y2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 1; % provare i valori 1 e 0.5
tspan = [0, 20];
odefun = @(t, y) [mu*y(1) - y(2) - y(1)*( y(1).^2 + 1.5 * y(2).^2 );
    y(1) + mu*y(2) - y(2)*(y(1).^2 + 0.5 * y(2).^2 )];

y0 = [3; -3];

[t, y] = ode15s(odefun, tspan, y0);

y1 = y(:, 1);
y2 = y(:, 2);

subplot(2, 2, 1)
plot(t, y1, LineWidth=2)
grid on
xlabel('tempo')
title('y1', Interpreter='latex')

subplot(2, 2, 2)
plot(t, y2, LineWidth=2)
grid on
xlabel('tempo')
title('y2', Interpreter='latex')

subplot(2, 2, 3)
r = sqrt(y1.^2 + y2.^2);
plot(t, r, LineWidth=2)
grid on
xlabel('tempo')
title('r', Interpreter='latex')

subplot(2, 2, 4)
plot(y1, y2, LineWidth=2)
grid on
title('Spazio delle fasi', Interpreter='latex')
xlabel('y1')
ylabel('y2')


%% punto b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fare plot delle funzioni u, u1, u2 in funzione del tempo e fare plot
% nello spazio delle fasi (u1, u2).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = [0, 100]; % provare i valori [0, 100] e [0, 300];

a = pi/3.1; % provare i valori pi/3.1 e sqrt(2)  

u = @(t) exp(-t) + sin(t) + cos(sqrt(2)*t);
u1 = @(t) exp(-t) + sin(t);
u2 = @(t) exp(-t) + cos(a*t);

figure;
subplot(2, 2, 1)
fplot(u, tspan, LineWidth=2);
grid on
title(['u, con a =', num2str(a)], Interpreter="latex")

subplot(2, 2, 2)
fplot(u1, tspan, LineWidth=2);
grid on
title(['u1, con a =', num2str(a)], Interpreter="latex")

subplot(2, 2, 3)
fplot(u2, tspan, LineWidth=2);
grid on
title(['u2, con a =', num2str(a)], Interpreter="latex")

subplot(2, 2, 4)
fplot(u1, u2, tspan, LineWidth=2);
grid on
title(['spazio delle fasi, con a =', num2str(a)], Interpreter="latex")
xlabel('u1')
ylabel('u2')


%% punto c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attrattore di Lorenz:
% Plottare y1, y2, y3, r in funzione del tempo e nello spazio
% tridimensionale (y1, y2, y3).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = [0, 300];
sigma = 10; 
r = 28; 
b = 8/3;

y0 = [0; 10; 0];
odefun = @(t, y) [sigma*(y(2)-y(1)); r*y(1) - y(2) - y(1)*y(3); y(1)*y(2) - b*y(3)];

[t, y] = ode15s(odefun, tspan, y0);
y1 = y(:, 1);
y2 = y(:, 2);
y3 = y(:, 3);

figure; 
subplot(4, 1, 1)
plot(t, y1, LineWidth=2);
grid on
xlabel('tempo')
title('y1', Interpreter='latex')

subplot(4, 1, 2)
plot(t, y2, LineWidth=2);
grid on
xlabel('tempo')
title('y2', Interpreter='latex')

subplot(4, 1, 3)
plot(t, y3, LineWidth=2);
grid on
xlabel('tempo')
title('y3', Interpreter='latex')

subplot(4, 1, 4)
r = sqrt(y1.^2 + y2.^2);
plot(t, r, LineWidth=2);
grid on
xlabel('tempo')
title('r', Interpreter='latex')


figure; 

plot3(y(:, 1), y(:, 2), y(:, 3), LineWidth=1)
grid on
xlabel('y1')
ylabel('y2')
zlabel('y3')
title('Attrattore di Lorentz')

%% punto d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Moto di un corpo in orbita intorno ad un altro:
% Plottare y1, y2, r in funzione del tempo e fare plot nello spazio delle 
% fasi (y1, y2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 0.25;
alpha = pi/4;

R =@(v, w) ((v.^2 + w.^2).^(1.5)) / alpha.^2;

tspan = [0, 25]; % provare anche con [0, 100] e [0, 300]

y0 = [1-rho; 0; 0; alpha*sqrt( (1+rho) / (1-rho) )];

odefun = @(t, y) [y(3); y(4); -y(1)/ R(y(1), y(2)); -y(2)/ R(y(1), y(2)) ];

[t, y] = ode15s(odefun, tspan, y0);

y1 = y(:, 1);
y2 = y(:, 2);
y3 = y(:, 3);
y4 = y(:, 4);

figure; 
subplot(2, 2, 1)
plot(t, y1, LineWidth=2)
grid on
xlabel('tempo')
title('y1', Interpreter='latex')

subplot(2, 2, 2)
plot(t, y2, LineWidth=2)
grid on
xlabel('tempo')
title('y2', Interpreter='latex')

subplot(2, 2, 3)
r = sqrt(y1.^2 + y2.^2);
plot(t, r, LineWidth=2)
grid on
xlabel('tempo')
title('r', Interpreter='latex')

subplot(2, 2, 4)
plot(y1, y2, LineWidth=2)
grid on
xlabel('tempo')
title('Spazio delle fasi y1-y2', Interpreter='latex')
xlabel('y1')
ylabel('y2')




