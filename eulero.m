function [t, ya] = eulero(f, tspan, y0, h)
% metodo eulero per esercitazione 2
% f è un puntatore a funzione
% tspan è un vettore [tempo iniziale, tempo finale]
% y0 valore iniziale
% h è uno scalare: ampiezza passo

N = ceil((tspan(2)-tspan(1)) /  h); 
t = linspace(tspan(1), tspan(2), N+1)';
ya = zeros(2, N+1);
ya(:, 1) = y0;

for n=1:N
    K1n = f(t(n), ya(:, n));
    ya(:, n+1) = ya(:, n) + h*K1n;
end

ya = ya.';

end