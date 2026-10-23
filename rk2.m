function [t, ya] = rk2(f, tspan, y0, h, alpha)
% metodo runge kutta per esercitazione 2
    N = ceil((tspan(2)-tspan(1)) /  h); 
    h = (tspan(2) - tspan(1)) / N;
    t = linspace(tspan(1), tspan(2), N+1)';

    ya = zeros(2, N+1);
    ya(:, 1) = y0;

    for n=1:N
        K1n = f(t(n), ya(:, n));
        K2n = f(t(n)+alpha*h, ya(:, n) + alpha * h * K1n);
        
        ya(:, n+1) = ya(:, n) + h * ( (1-1/(2*alpha))*K1n + 1/(2*alpha) * K2n );

    end

    ya = ya.';
    
end