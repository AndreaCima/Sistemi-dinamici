function [t, ya] = rk4(f, tspan, y0, h)
% metodo Runge Kutta di ordine 4 per esercitazione 2

    if isscalar(h)
        N = ceil((tspan(2)-tspan(1)) /  h); 
        h = (tspan(2) - tspan(1)) / N * ones(N, 1).';
        t = linspace(tspan(1), tspan(2), N+1)';
    else % h è un vettore
        N = length(h); 
        h = h(:); % per essere sicuro che sia un vettore colonna
        t =  [tspan(1); cumsum(h)];
    end
    
    
    d = length(y0);
    ya = zeros(d, N+1);
    ya(:, 1) = y0(:);

    for n = 1:N
        K1n = f(t(n), ya(:, n));
        K2n = f(t(n) + h(n)/2, ya(:, n) + h(n)/2 * K1n);
        K3n = f(t(n) + h(n)/2, ya(:, n) + h(n)/2 * K2n);
        K4n = f(t(n) + h(n), ya(:, n) + h(n)*K3n);

        ya(:, n+1) = ya(:, n) + h(n)/3 * (K1n / 2 + K2n + K3n + K4n / 2);
    end

    ya = ya.';
end