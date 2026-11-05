% plot regioni di assoluta stabilità
clear; clc; close all

%% Eulero esplicito
x_range= linspace(-3, 1);
y_range = linspace(-2, 2);

[X, Y] = meshgrid(x_range, y_range);

Z = X + 1i*Y; % punti nel piano complesso

in = Z( abs(Z + 1) < 1 );
out = Z( abs(Z + 1) >= 1 );

figure;
scatter(real(in), imag(in), 'ob', DisplayName = 'in')
hold on
scatter(real(out), imag(out), 'xr', DisplayName = 'out')
grid on
axis square
title('Eulero esplicito', Interpreter='latex')
legend(Location="best")

figure; 
abs_Z = min(abs(Z + 1), 1);
surf(X, Y, abs_Z) 
colorbar
axis square
title('Eulero esplicito', Interpreter='latex')

figure; 
contour(X, Y, abs_Z)
grid on
axis square
title('Eulero esplicito', Interpreter='latex')

%% Eulero implicito

x_range= linspace(-1, 3);
y_range = linspace(-2, 2);

[X, Y] = meshgrid(x_range, y_range);

Z = X + 1i*Y; % punti nel piano complesso

in = Z( abs(Z - 1) > 1 );
out = Z( abs(Z - 1) <= 1 );

figure;
scatter(real(in), imag(in), 'ob', DisplayName = 'in')
hold on
scatter(real(out), imag(out), 'xr', DisplayName = 'out')
grid on
axis square
title('Eulero implicito', Interpreter='latex')

figure; 
abs_Z = min(1./abs(Z - 1), 1);
surf(X, Y, abs_Z) 
colorbar
axis square
title('Eulero implicito', Interpreter='latex')

figure; 
contour(X, Y, abs_Z)
grid on
axis square
title('Eulero implicito', Interpreter='latex')

%% Runge Kutta 2

A_rk2 = [0 0; 1 0];
b_rk2 = [0.5, 0.5];
c_rk2 = [0; 1];

x_range= linspace(-3, 1);
y_range = linspace(-2, 2);

[X, Y] = meshgrid(x_range, y_range);

Z = X + 1i*Y; % punti nel piano complesso

rZ = r(A_rk2, b_rk2, Z);
in = Z( abs(rZ) < 1 );
out = Z( abs(rZ) >= 1 );

figure;
scatter(real(in), imag(in), 'ob', DisplayName = 'in')
hold on
scatter(real(out), imag(out), 'xr', DisplayName = 'out')
grid on
axis square
title('RK 2', Interpreter='latex')

figure; 
abs_Z = min(abs(rZ), 1);
surf(X, Y, abs_Z) 
colorbar
axis square
title('RK 2', Interpreter='latex')

figure; 
contour(X, Y, abs_Z)
grid on
axis square
title('RK 2', Interpreter='latex')

%% Runge Kutta 4

A_rk4 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
b_rk4 = [1/6 1/3 1/3 1/6];
c_rk4 = [0 1/2 1/2 1];

x_range= linspace(-3, 1);
y_range = linspace(-3, 3);

[X, Y] = meshgrid(x_range, y_range);

Z = X + 1i*Y; % punti nel piano complesso

rZ = r(A_rk2, b_rk2, Z);

in = Z( abs(rZ) < 1 );
out = Z( abs(rZ) >= 1 );

figure;
scatter(real(in), imag(in), 'ob', DisplayName = 'in')
hold on
scatter(real(out), imag(out), 'xr', DisplayName = 'out')
grid on
axis square
title('RK 4', Interpreter='latex')

figure; 
abs_Z = min(abs(rZ), 1);
surf(X, Y, abs_Z) 
colorbar
axis square
title('RK 4', Interpreter='latex')

figure; 
contour(X, Y, abs_Z)
grid on
axis square
title('RK 4', Interpreter='latex')


function rZ = r(A, b, Z)
    rZ = zeros(size(Z));
    dim = length(b);
    for h = 1:size(Z, 1)
        for k = 1:size(Z, 2)
            z = Z(h, k);
            rz = 1 + z*(b*( (eye(dim)-z*A) \ ones(dim, 1) ));
            rZ(h, k) = rz;
        end
    end
end
