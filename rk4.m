function [t, ya] = rk4(f, tspan, y0, h)
% metodo Runge Kutta di ordine 4 per esercitazione 2
    
    N = ceil((tspan(2)-tspan(1)) /  h); 
    t = linspace(tspan(1), tspan(2), N+1)';
    
    ya = zeros(2, N+1);
    ya(:, 1) = y0;

    for n = 1:N
        K1n = f(t(n), ya(:, n));
        K2n = f(t(n) + h/2, ya(:, n) + h/2 * K1n);
        K3n = f(t(n) + h/2, ya(:, n) + h/2 * K2n);
        K4n = f(t(n) + h, ya(:, n) + h*K3n);

        ya(:, n+1) = ya(:, n) + h/3 * (K1n/2 + K2n + K3n + K4n/2);
    end

    ya = ya.';
end