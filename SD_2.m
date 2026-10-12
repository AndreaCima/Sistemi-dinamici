clear; clc; close all; 

y01 = 1; 
y02 = -3;
mu = 1;
tspan = [0, 20];
alpha = 1; % provare 1 e 0.5

y0 = [y01; y02];

f = @(t, y) [ ...
    mu*y(1) + y(2) - (y(1) + y(2))*(y(1)^2 + y(2)^2) ; ...  
    -y(1) + mu*y(2) - (y(2) - y(1))*(y(1)^2 + y(2)^2)   ...  
    ];

h = 0.001;

% calcolo approssimazioni con metodi di eulero e runge-kutta
[t, ya_eulero] = eulero(f, tspan, y0, h);
[~, ya_rk2] = rk2(f, tspan, y0, h, alpha);

ye = esatta(t, mu, y0); % soluzione esatta

% errori

er1_eulero = ye(:, 1) - ya_eulero(:, 1);
er2_eulero = ye(:, 2) - ya_eulero(:, 2);
er1_rk2 = ye(:, 1) - ya_rk2(:, 1);
er2_rk2 = ye(:, 2) - ya_rk2(:, 2);

Einf_eulero = max(abs(ye-ya_eulero), [], 2);
Einf_rk2 = max(abs(ye-ya_rk2), [], 2);

E2_eulero = sqrt(sum( (ye-ya_eulero).^2, 2));
E2_rk2 = sqrt(sum( (ye-ya_rk2).^2, 2));

nye = sqrt(sum(ye.^2, 2));
nya_eulero = sqrt(sum(ya_eulero.^2, 2));
nya_rk2 = sqrt(sum(ya_rk2.^2, 2));

err_eulero = norm(E2_eulero)/norm(nye);
err_rk2 = norm(E2_rk2)/norm(nye);


figure; 
subplot(2, 2, 1)
plot(t, ya_eulero(:, 1), 'r', DisplayName='y1 eulero')
hold on
plot(t, ye(:, 1), '--b', DisplayName='y1 esatta')
plot(t, ya_eulero(:, 2), 'r', DisplayName='y2 eulero')
plot(t, ye(:, 2), '--b', DisplayName='y2 esatta')
legend(Location="best")
grid on
xlabel('tempo')


subplot(2, 2, 2)
plot(t, er1_eulero, DisplayName='er1 eulero')
hold on
plot(t, er2_eulero, DisplayName='er2 eulero')
plot(t, Einf_eulero, DisplayName='Einf eulero')
plot(t, E2_eulero, DisplayName='E2 eulero')
legend(Location="best")
grid on
xlabel('tempo')


subplot(2, 2, 3)

plot(t, nya_eulero, 'r', DisplayName='errore nya eulero')
hold on
plot(t, nye, '--b', DisplayName='errore nye')
legend(Location="best")
grid on
xlabel('tempo')

subplot(2, 2, 4)
plot(ya_eulero(:, 1), ya_eulero(:, 2), 'r', DisplayName='approssimazione eulero')
hold on
plot(ye(:, 1), ye(:, 2), '--b',  DisplayName='soluzione esatta')
legend(Location="best")
axis equal
grid on
title('spazio delle fasi', Interpreter='latex')


figure; 
subplot(2, 2, 1)
plot(t, ya_rk2(:, 1), 'r', DisplayName='y1 rk2')
hold on
plot(t, ye(:, 1), '--b', DisplayName='y1 esatta')
plot(t, ya_rk2(:, 2), 'r', DisplayName='y2 rk2')
plot(t, ye(:, 2), '--b', DisplayName='y2 esatta')
legend(Location="best")
xlabel('tempo')
grid on


subplot(2, 2, 2)
plot(t, er1_rk2, DisplayName='er1 rk2')
hold on
plot(t, er2_rk2, DisplayName='er2 rk2')
plot(t, Einf_rk2, DisplayName='Einf v')
plot(t, E2_rk2, DisplayName='E2 rk2')
legend(Location="best")
grid on
xlabel('tempo')


subplot(2, 2, 3)
plot(t, nya_rk2, 'r', DisplayName='errore nya rk2')
hold on
plot(t, nye, '--b', DisplayName='errore nye')
legend(Location="best")
grid on
xlabel('tempo')

subplot(2, 2, 4)
plot(ya_rk2(:, 1), ya_rk2(:, 2), 'r', DisplayName='approssimazione rk2')
hold on
plot(ye(:, 1), ye(:, 2), '--b', DisplayName='soluzione esatta')
legend(Location="best")
axis equal
grid on
title('spazio delle fasi', Interpreter='latex')



function ye = esatta(t, mu, y0)
    y01 = y0(1); 
    y02 = y0(2);
    
    r0 = sqrt(y01^2 + y02^2);
    theta0 = atan2(y02, y01);

    phi = @(t) r0^2 + (mu - r0^2) * exp(-2*mu*t);
    r = @(t) (sqrt(mu)*r0) ./ sqrt(phi(t));
    theta = @(t) theta0 - t + mu * (t + (log(phi(t)) - log(mu)) / (2*mu));
    
    % soluzioni esatte
    y1=@(t) r(t) .* cos(theta(t));
    y2=@(t) r(t) .* sin(theta(t)); 

    assert(norm(y0 - [y1(t(1)); y2(t(1))]) < 1e-5, 'I dati iniziali non coincidono')

    ye = [y1(t) y2(t)];
end

