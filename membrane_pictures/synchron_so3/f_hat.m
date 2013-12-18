function f = f_hat(l, image1)

[n1, n2] = size(image1);

scale = 4;
theta_vec = linspace(pi/scale, pi - pi/scale, n1);
phi_vec = linspace(pi/scale, pi - pi/scale, n2);

theta = repmat(theta_vec', 1, n2);
phi = repmat(phi_vec, n1, 1);

size_vec = size(theta);
theta = reshape(theta, 1, []);
phi = reshape(phi, 1, []);
image1 = reshape(image1, 1, []);

f = (Y_lm(l, theta, phi)*image1').';
% f = sum(Y_lm(l, theta, phi).*repmat(image1, 2*l+1, 1), 2)'
% y = Y_lm(l, theta, phi);
% f = zeros(1, 2*l+1);
% for i=1:2*l+1
%     f(i) = sum(y(i,:).*image1);
% end
