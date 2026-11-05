clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y0 = [0.1; 1];
mu = -5;
T = 10; 
h = 1e-1;

A = [mu 0; 1+mu -1];
tspan = [0 T];
f = @(t, y) A*y;

% soluzione esatta ricavata analiticamente
esatta = @(t) [0.1*exp(mu*t) 0.9*exp(-t)+0.1*exp(mu*t)];

[t, ya] = rk4(f, tspan, y0, h);
ye = esatta(t);

er1 = abs(ye(:, 1)-ya(:, 1));
er2 = abs(ye(:, 2)-ya(:, 2));

ye_norm2 = sqrt(sum(ye.^2, 2));
ya_norm2 = sqrt(sum(ya.^2, 2));

tiledlayout(2, 2)
title(tiledlayout(2, 2), sprintf('h = %.2f', h), Interpreter='latex')

nexttile
plot(t, ya(:, 1), 'b', LineWidth=2, DisplayName='appr1')
hold on
grid on
plot(t, ya(:, 2), 'r', LineWidth=2, DisplayName='appr2')
plot(t, ye(:, 1), '--o', LineWidth=2, DisplayName='ex1')
plot(t, ye(:, 2), '--o', LineWidth=2, DisplayName='ex2')
legend(Location="best")

nexttile 
plot(t, er1, LineWidth=2, DisplayName='er1')
hold on
grid on
plot(t, er2, '--', LineWidth=2, DisplayName='er2')
legend(Location="best")
xlabel('tempo')


nexttile 
plot(t, ye_norm2, LineWidth=2, DisplayName='$\|ya\|_2$')
hold on
grid on
plot(t, ya_norm2, '--', LineWidth=2, DisplayName='$\|ya\|_2$')
xlabel('tempo')
legend(Location="best", Interpreter="latex")

nexttile
plot(ya(:, 1), ya(:, 2), LineWidth=2, DisplayName='soluzione approssimata')
grid on
hold on
plot(ye(:, 1), ye(:, 2), '--', LineWidth=2, DisplayName='soluzione esatta')
title('Spazio delle fasi', Interpreter='latex')
xlabel('prima componente')
ylabel('seconda componente')
legend(Location="best")




%% passo h variabile
clear ya ye
h_test = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];

treshold = 4/abs(mu);

for idx = 1:length(h_test)
    h = h_test(idx);
    N1 = ceil(treshold/0.1);
    N2 = ceil((T- treshold)/h);
    h_vec = [0.1*ones(N1, 1); h*ones(N2, 1)];

    [t, ya] = rk4(f, tspan, y0, h_vec);
    ye = esatta(t);

    er1 = abs(ye(:, 1)-ya(:, 1));
    er2 = abs(ye(:, 2)-ya(:, 2));

    ye_norm2 = sqrt(sum(ye.^2, 2));
    ya_norm2 = sqrt(sum(ya.^2, 2));

    figure;
    tiledlayout(2, 2)
    title(tiledlayout(2, 2), sprintf('h1 = 0.1, h2 = %.2f', h), Interpreter='latex')
    
    nexttile
    plot(t, ya(:, 1), 'b', LineWidth=2, DisplayName='appr1')
    hold on
    grid on
    plot(t, ya(:, 2), 'r', LineWidth=2, DisplayName='appr2')
    plot(t, ye(:, 1), '--o', LineWidth=2, DisplayName='ex1')
    plot(t, ye(:, 2), '--o', LineWidth=2, DisplayName='ex2')
    legend(Location="best")
    
    nexttile 
    plot(t, er1, LineWidth=2, DisplayName='er1')
    hold on
    grid on
    plot(t, er2, '--', LineWidth=2, DisplayName='er2')
    legend(Location="best")
    xlabel('tempo')
    
    
    nexttile 
    plot(t, ye_norm2, LineWidth=2, DisplayName='$\|ya\|_2$')
    hold on
    grid on
    plot(t, ya_norm2, '--', LineWidth=2, DisplayName='$\|ya\|_2$')
    xlabel('tempo')
    legend(Location="best", Interpreter="latex")
    
    nexttile
    plot(ya(:, 1), ya(:, 2), LineWidth=2, DisplayName='soluzione approssimata')
    grid on
    hold on
    plot(ye(:, 1), ye(:, 2), '--', LineWidth=2, DisplayName='soluzione esatta')
    title('Spazio delle fasi', Interpreter='latex')
    xlabel('prima componente')
    ylabel('seconda componente')
    legend(Location="best")


end

% dai diversi plot si vede che il valore massimo di h tale per cui si ha
% ancora una buona approssimazione è h = 0.5

Cols = {'h', 'er1', 'er2', 'ye2', 'ya2'};
Cols_type = {'double', 'double', 'double', 'double', 'double'};

Table = table('Size', [length(h_test), length(Cols)], ...
    'VariableTypes', Cols_type, ...
    'VariableNames', Cols);

N1 = ceil(treshold/0.1);
N2 = ceil((T- treshold)/0.5); % uso il valore h = 0.5
N = N1 + N2;
h_vec = [0.1*ones(N1, 1); 0.5*ones(N2, 1)];

[t, ya] = rk4(f, tspan, y0, h_vec);
ye = esatta(t);

er1 = abs(ye(:, 1)-ya(:, 1));
er2 = abs(ye(:, 2)-ya(:, 2));

ye_norm2 = sqrt(sum(ye.^2, 2));
ya_norm2 = sqrt(sum(ya.^2, 2));

figure;
tiledlayout(2, 2)
title(tiledlayout(2, 2), 'h1 = 0.1, h2 = 0.5', Interpreter='latex')

nexttile
plot(t, ya(:, 1), 'b', LineWidth=2, DisplayName='appr1')
hold on
grid on
plot(t, ya(:, 2), 'r', LineWidth=2, DisplayName='appr2')
plot(t, ye(:, 1), '--o', LineWidth=2, DisplayName='ex1')
plot(t, ye(:, 2), '--o', LineWidth=2, DisplayName='ex2')
legend(Location="best")

nexttile
plot(t, er1, LineWidth=2, DisplayName='er1')
hold on
grid on
plot(t, er2, '--', LineWidth=2, DisplayName='er2')
legend(Location="best")
xlabel('tempo')


nexttile
plot(t, ye_norm2, LineWidth=2, DisplayName='$\|ya\|_2$')
hold on
grid on
plot(t, ya_norm2, '--', LineWidth=2, DisplayName='$\|ya\|_2$')
xlabel('tempo')
legend(Location="best", Interpreter="latex")

nexttile
plot(ya(:, 1), ya(:, 2), LineWidth=2, DisplayName='soluzione approssimata')
grid on
hold on
plot(ye(:, 1), ye(:, 2), '--', LineWidth=2, DisplayName='soluzione esatta')
title('Spazio delle fasi', Interpreter='latex')
xlabel('prima componente')
ylabel('seconda componente')
legend(Location="best")



