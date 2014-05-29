function y = wigner_semicircle(x, R, alpha)

y = alpha * sqrt(R^2-x.^2);

y(abs(x) > R) = 0;