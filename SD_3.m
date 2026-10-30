clear; clc; close all

y0 = [0.1; 1];
mu = -5;
T = 10; 
h_cost = 1e-3;

A = [mu 0; 1+mu -1];
tspan = [0 T];
f = @(t, y) A*y;


esatta =@(t) [0.1*exp(mu*t) 0.9*exp(-t)+0.1*exp(mu*t)];

[t, ya_h_cost] = rk4(f, tspan, y0, h_cost);
ye = esatta(t);


plot(t, ye, '-b', DisplayName='soluzione esatta')
hold on
grid on
plot(t, ya_h_cost, '--r', DisplayName='approssimazione')
xlabel('tempo', Interpreter='latex')
legend(Location='best')
title('passo h costante')


% passo h variabile 

h_test = [0.1 0.2 0.3 0.4 0.5];




for h = h_test
   
    h_vec = h_fun(h, t, mu);
    [t, ya_h_var] = rk4(f, tspan, y0, h_vec);
    
    figure;
    ye = esatta(t);
    plot(t, ye, '-b', DisplayName='soluzione esatta')
    
    grid on
    hold on 
    
    plot(t, ya_h_var, '--', DisplayName='approssimazione')
    title(sprintf('h = %g', h), Interpreter="latex")
    legend(Location="best")

end

