function y = Y(l, theta, phi)

if size(theta) > 2 || size(theta, 1) ~= 1 || size(phi) > 2 || size(phi, 1) ~= 1
    disp('ERROR: inconsistent dimension')
end

P = legendre(l, cos(theta));

m = repmat((-l:1:-1)', 1, size(P,2));
P = [(-1).^(-m).*factorial(l+m)./factorial(l-m).*P(end:-1:2,:); P];

m = repmat((-l:l)', 1, size(P,2));
y = sqrt(((2*l+1)*factorial(l-m))./(4*pi*factorial(l+m))) .* P .* exp(sqrt(-1)*m.*repmat(phi, 2*l+1, 1));