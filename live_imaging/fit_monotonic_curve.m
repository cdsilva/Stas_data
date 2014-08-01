function zopt = fit_monotonic_curve(y, slack)

n = length(y);

fun = @(z) sum((z-y).^2);

z0 = linspace(y(1), y(end), n)';

A = zeros(n-1, n);
for i=1:n-1
    A(i, i) = 1;
    A(i, i+1) = -1;
end
b = slack * ones(n-1, 1);

options = optimoptions(@fmincon, 'display','off', 'algorithm', 'interior-point');
zopt = fmincon(fun, z0, A, b, [], [], [], [], [], options);



