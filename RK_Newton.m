function [t, ya, it] = RK_Newton(A, b, c, f, df, y0, T, nsteps)
    % calcola approssimazione nuemrica di una ODE con Runge Kutta implicito
    % o semi-implicito usando metodo di Newton ad ognni iterazione.

    % parametri del metodo di Newton 
    maxiter = 20; 
    tol = 1e-10;

    dispProgress = false;

    t = linspace(0, T, nsteps)';
    h = t(2)-t(1);
    ya = zeros(nsteps, 1);
    ya(1) = y0; 
    nstages = length(b);

    it = zeros(nsteps-1, 1);
    ERR = zeros(nsteps-1, 1);

    for n = 1:nsteps - 1
        if dispProgress
            fprintf('Progresso %i / %i \t', n, nsteps-1)
        end
        tn = t(n);
        yn = ya(n);

        K = zeros(nstages, 1);

        iter = 0; 
        err = 1;

        while (err> tol) && (iter <maxiter)
            iter = iter +1; 

            Y_stage = yn + h*(A*K);
            T_stage = tn + h*c;

            f_val = f(T_stage, Y_stage);
            df_val = df(T_stage, Y_stage);

            Phi = K - f_val;
            DF = eye(nstages) - h*(diag(df_val) * A);

            dK = - (DF\Phi);

            K = K + dK;
            err= max(abs(dK));
        end
        ya(n+1) = yn + h*(b*K);
        it(n) = iter; % salvo il numero di iterazioni
        if dispProgress
            fprintf('iterazioni richieste = %i \n', iter)
        end
        ERR(n) = err;

    end


end