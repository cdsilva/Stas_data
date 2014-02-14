function im1 = rotate_image(R, im2, angle_proj)

npixels = size(im2, 1);

[theta, dxpixels, dypixels] = find_specs(R, npixels, angle_proj);

im2 = imrotate(im2, theta, 'crop');
im2 = circshift(im2, [dxpixels dypixels]);

im1 = im2;

function [theta, dx, dy] = find_specs(R, npixels, angle_proj)

alpha = atan2(R(3,2), R(3,3));
gamma = atan2(R(2,1), R(1,1));
beta = -asin(R(3,1));
%beta = atan2(-R(3,1), R(2,1)/sin(gamma));

theta = alpha * 180 / pi;
dx = round(beta * npixels / angle_proj);
dy = round(gamma * npixels / angle_proj);


