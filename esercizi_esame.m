%% 19/02/2018 es 4

EPSILON = [1e-2 1e-4 1e-6];
y_10 = 1.5; y_20 = 1;
T = 10; 

tspan = [0 T];

y0 = [y_10; y_20];



Cols = {'metodo', 'epsilon', 'n_step', 'err'}; 
Cols_type = {'string', 'double', 'double', 'double'};
Table = table('Size', [length(EPSILON)*3, length(Cols)], ...
    'VariableTypes', Cols_type, ...
    'VariableNames', Cols);



for idx = 1:length(EPSILON)
    
    epsilon = EPSILON(idx);
    fprintf('epsilon = %f --> ', epsilon)

    f = @(t, y) [1/epsilon * (-y(1) + y(2)./(1+y(2)) ) - y(2) ./(1+y(2)).^2; ...
        -y(2)];

    ye1 = @(t) (y_10 - y_20 ./(1+y_20) ).*exp(-t/epsilon) + (y_20 *exp(-t)) ./ (1+y_20*exp(-t));
    ye2 = @(t) y_20 * exp(-t);
    
    fprintf('ode15s \t')
    [t, ya] = ode15s(f, tspan, y0);
    ye = [ye1(t) ye2(t)];

    E2 = sqrt(sum((ye-ya).^2,2));
    nye = sqrt(sum(ye.^2,2));
    err = norm(E2)/norm(nye);

    Table.metodo(3*(idx-1)+1) = 'ode15s';
    Table.epsilon(3*(idx-1)+1) = epsilon;
    Table.n_step(3*(idx-1)+1) = length(t);
    Table.err(3*(idx-1)+1) = err;

    fprintf('ode23 \t')
    [t, ya] = ode23(f, tspan, y0);
    ye = [ye1(t) ye2(t)];

    E2 = sqrt(sum((ye-ya).^2,2));
    nye = sqrt(sum(ye.^2,2));
    err = norm(E2)/norm(nye);

    Table.metodo(3*(idx-1)+2) = 'ode23';
    Table.epsilon(3*(idx-1)+2) = epsilon;
    Table.n_step(3*(idx-1)+2) = length(t);
    Table.err(3*(idx-1)+2) = err;

    fprintf('ode45 \n')
    [t, ya] = ode45(f, tspan, y0);
    ye = [ye1(t) ye2(t)];

    E2 = sqrt(sum((ye-ya).^2,2));
    nye = sqrt(sum(ye.^2,2));
    err = norm(E2)/norm(nye);

    Table.metodo(3*(idx-1)+3) = 'ode45';
    Table.epsilon(3*(idx-1)+3) = epsilon;
    Table.n_step(3*(idx-1)+3) = length(t);
    Table.err(3*(idx-1)+3) = err;

end

disp(Table)


%% 20/06/2018 es 2

Phi = @(y) y(:, 1).^2 + y(:, 2).^2 + y(:, 3).^2 -1;
gradPhi = @(y) [2*y(1); 2*y(2); 2*y(3)];
HPhi = 2*eye(3);
F = [0 0 -9.8]';

lambda = @(y) (2*(y(4).^2 + y(5).^2 + y(6).^2) + gradPhi(y).' * F) ./ (gradPhi(y).'*gradPhi(y));

f = @(t, y) [y(4); y(5); y(6); ...
    F(1) - 2*lambda(y) * y(1); ...
    F(2) - 2*lambda(y) * y(2); ...
    F(3) - 2*lambda(y) * y(3)];

T = 25; 
y0 = [0 1 0 0.8 0 1.2]';

tspan = [0, T];

options = odeset('RelTol', 1e-5);

[t, ya] = ode15s(f, tspan, y0, options);

ya1 = ya(:, 1);
ya2 = ya(:, 2);
ya3 = ya(:, 3);
ya4 = ya(:, 4);
ya5 = ya(:, 5);
ya6 = ya(:, 6);

figure; 
plot3(ya1, ya2, ya3)
grid on
hold on
xlabel('ya1')
ylabel('ya2')
zlabel('ya3')
% scatter3(y0(1), y0(2), y0(3))
axis equal

% plot sfera
[X, Y, Z] = sphere(50);
surf(X, Y, Z, FaceAlpha=0.1, EdgeColor='none')


