function [t, ya] = butcher_esplicito(A, b, c, h, tspan, f, y0)
% Calcolo approssimazione ODE con metodi Runge Kutta espliciti, basandosi
% sulla tabella di Butcher
arguments (Input)
    A % matrice s x s
    b % vettore 1 x s
    c % vettore s x 1
    h % scalare, passo temporale h_{n+1} - h_n
    tspan % intervallo temporale, 2 x 1
    f % funzione, dato del problema
    y0 % vettore, dato iniziale del problema

end

arguments (Output)
    t % vettore dei tempi
    ya % soluzione approssimata
end

N = ceil((tspan(2)-tspan(1)) /  h); 
h = (tspan(2) - tspan(1)) / N;
t = linspace(tspan(1), tspan(2), N+1)';
s = length(c); % numero stadi

ya = zeros(length(y0), N+1);
ya(:, 1) = y0;


for n = 1:N
    T = t(n) + c*h; % istanti temporali t_{si}
    
    K = zeros(s, length(y0));
    for i = 1:s
        K(i, :) = f(T(i), ya(:, n) + h * ( A(i, :)*K ).');
    end
    
    ya(:, n+1) = ya(:, n) + h * ((b*K).'); % approssimazione al passo t_{n+1}

end

ya = ya.';

end