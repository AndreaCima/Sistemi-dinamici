clear; clc; close all;

f = @(t, y) y - t + 1;
df = @(t, y) ones(size(y)); % derivata parziale rispetto ad y

ye = @(t) exp(t) + t; % soluzione esatta, ricavata in modo analitico

y0 = 1; % condizione iniziale
T = 1; % studio il problema per t \in [0, T]
nsteps = 41; % numero passi temporali

% criteri di arresto per metodo di Newton
maxiter = 20; 
tol = 1e-10; 

%%%%%%%%%%%%%%%%%%%%%% Tabelle di Butcher %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eulero implicito
A_EI = 1; 
b_EI = 1; 
c_EI = 1; 

% trapezi
A_tr = [0 0; 0.5 0.5];
b_tr = [0.5 0.5];
c_tr = [0; 1];

% % RK 2
A_RK2 = [0 0; 1 0];
b_RK2 = [0.5, 0.5];
c_RK2 = [0; 1];


[t, ya, ~] = RK_Newton(A_EI, b_EI, c_EI, f, df, y0, T, nsteps);

figure; 
plot(t, ye(t), LineWidth=2, DisplayName='soluzione esatta')
grid on
hold on
plot(t, ya, LineWidth=2, DisplayName='approssimazione')
legend(Location="best", Interpreter="latex")
xlabel('tempo', Interpreter='latex')
title('Eulero implicito', Interpreter='latex')


[t, ya, ~] = RK_Newton(A_tr, b_tr, c_tr, f, df, y0, T, nsteps);

figure; 
plot(t, ye(t), LineWidth=2, DisplayName='soluzione esatta')
grid on
hold on
plot(t, ya, LineWidth=2, DisplayName='approssimazione')
legend(Location="best", Interpreter="latex")
xlabel('tempo', Interpreter='latex')
title('Trapezi', Interpreter='latex')


[t, ya, ~] = RK_Newton(A_RK2, b_RK2, c_RK2, f, df, y0, T, nsteps);

figure; 
plot(t, ye(t), LineWidth=2, DisplayName='soluzione esatta')
grid on
hold on
plot(t, ya, LineWidth=2, DisplayName='approssimazione')
legend(Location="best", Interpreter="latex")
xlabel('tempo', Interpreter='latex')
title('Runge-Kutta 2', Interpreter='latex')

%% Studio ordine di convergenza

NSTEPS = 2.^(4:12)';

ERR = zeros(length(NSTEPS), 3);
ITER = zeros(length(NSTEPS), 3);
H = zeros(length(NSTEPS), 3);
for idx = 1:length(NSTEPS)
    nsteps = NSTEPS(idx);
    
    % Eulero implicito
    [t, ya, it] = RK_Newton(A_EI, b_EI, c_EI, f, df, y0, T, nsteps);
    ERR(idx, 1) = norm(ye(t)-ya, 'inf');
    ITER(idx, 1) = mean(it);
    H(idx, 1) = t(2)-t(1);
    
    % Trapezi
    [t, ya, it] = RK_Newton(A_tr, b_tr, c_tr, f, df, y0, T, nsteps);
    ERR(idx, 2) = norm(ye(t)-ya, 'inf');
    ITER(idx, 2) = mean(it);
    H(idx, 2) = t(2)-t(1);
    
    % Runge-Kutta 2
    [t, ya, it] = RK_Newton(A_RK2, b_RK2, c_RK2, f, df, y0, T, nsteps);
    ERR(idx, 3) = norm(ye(t)-ya, 'inf');
    ITER(idx, 3) = mean(it);
    H(idx, 3) = t(2)-t(1);
    

end


figure; 
PLOT = loglog(H, ERR, LineWidth=2);
PLOT(1).DisplayName = "Eulero implicito";
PLOT(1).Color = 'b';
PLOT(2).DisplayName = "Trapezi";
PLOT(2).Color = 'r';
PLOT(3).DisplayName = "RK2";
PLOT(3).Color = 'm';

grid on
hold on

loglog(H(:, 1), ERR(end, 1)*( H(:, 1)./H(end, 1) ).^1, '--ob', DisplayName='$O(h)$')
loglog(H(:, 2), ERR(end, 2)*( H(:, 2)./H(end, 2) ).^2, '--or', DisplayName='$O(h^2)$')
loglog(H(:, 3), ERR(end, 3)*( H(:, 3)./H(end, 3) ).^2, '--om', DisplayName='$O(h^2)$')

legend(Location="best", Interpreter="latex")

title('Ordini di convergenza con numero medio di iterazioni', Interpreter='latex')

for idx = 1:length(H)
    text(H(idx, 1), ERR(idx, 1), sprintf("%i", ITER(idx, 1)), ...
        VerticalAlignment='top', ...
        HorizontalAlignment='center')

    text(H(idx, 2), ERR(idx, 2), sprintf("%i", ITER(idx, 2)), ...
        VerticalAlignment='top', ...
        HorizontalAlignment='center')

    text(H(idx, 3), ERR(idx, 3), sprintf("%i", ITER(idx, 3)), ...
        VerticalAlignment='bottom', ...
        HorizontalAlignment='center')
end


