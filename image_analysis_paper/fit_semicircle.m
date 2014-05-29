function [alpha, R] = fit_semicircle(x, y)

sse = @(z) sum((y - wigner_semicircle(x, z(2), z(1))).^2);

fopt = fminsearch(sse, [1, 1]);

alpha = fopt(1);
R = fopt(2);
