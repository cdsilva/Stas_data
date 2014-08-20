function im1 = rotate_image(im2, R)

npixels = size(im2, 1);
angle_proj = pi/8;

dim = size(R, 1);

if dim == 3
    [theta, dx, dy] = find_specs_3d(R, angle_proj);
else % dim == 2
    theta = find_specs_2d(R);
    dx = 0;
    dy = 0;
end

RI = imref2d([npixels npixels],[-1 1],[-1 1]);

A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; dx dy 1]);
im1 = imwarp(im2, RI, A, 'outputview', RI);

function [theta, dx, dy] = find_specs_3d(R, angle_proj)

alpha = atan2(R(3,2), R(3,3));
gamma = atan2(R(2,1), R(1,1));
beta = -asin(R(3,1));
%beta = atan2(-R(3,1), R(2,1)/sin(gamma));

theta = alpha * 180 / pi;
dx = round(beta / angle_proj);
dy = round(gamma / angle_proj);

function theta = find_specs_2d(R)

alpha = atan2(R(2,1), R(2,2));
theta = alpha * 180 / pi;


