clear; clc; close all;

f = @(t, y) y - t + 1;
df = @(t, y) 1;

y0 = 1; 
T = 1; 
nsteps = 40; 
maxiter = 20; 
tol = 1e-10; 

%%%%%%%%%%%%%%%%%%%%%% tabelle di Butcher %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eulero implicito
A = 1; 
b = 1; 
c = 1; 

% trapezi
% A = [0 0; 0.5 0.5];
% b = [0.5 0.5];
% c = [0; 1];

% % RK 2
% A = [0 0; 1 0];
% b = [0.5, 0.5];
% c = [0; 1];



% altri parametri

t = linspace(0, T, nsteps);
h = t(2)-t(1);
Y = zeros(nsteps, 1);
Y(1) = y0; 
nstages = length(b); 

for n = 1:nsteps -1
    Kold = zeros(nstages, 1);
    iter = 0; 
    err = 1;

    while (err> tol) && (iter <maxiter)
        z = Y(n);
        iter = iter + 1;
    end

end

