function [x, y, z, f] = project_full_image(im1)


scale = 4;

[nim1, nim2] = size(im1);
im1 = double(imresize(uint8(im1), floor([nim1 nim2]/(scale-2))*(scale-2)));
[nim1, nim2] = size(im1);

%npoints = 20;
npoints1 = nim1/(scale-2);
npoints2a = nim2/(scale-2);
npoints2b = nim2*(scale+1)/(scale-2);

theta_vec = linspace(pi/scale, pi - pi/scale, nim1);
theta_vec1 = linspace(0, pi/scale, npoints1+1);
theta_vec1 = theta_vec1(1:end-1);
theta_vec2 = linspace(pi - pi/scale, pi, npoints1+1);
theta_vec2 = theta_vec2(2:end);

phi_vec = linspace(pi/scale, pi - pi/scale, nim2);
phi_vec1 = linspace(0, pi/scale, npoints2a+1);
phi_vec1 = phi_vec1(1:end-1);
phi_vec2 = linspace(pi - pi/scale, 2*pi, npoints2b+1);
phi_vec2 = phi_vec2(2:end);

theta = repmat([theta_vec1 theta_vec theta_vec2]', 1, nim2+npoints2a+npoints2b);
phi = repmat([phi_vec1 phi_vec phi_vec2], nim1+2*npoints1, 1);
f = zeros(size(theta));
f((npoints1+1):(npoints1+nim1),(npoints2a+1):(npoints2a+nim2)) = im1;

x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);

