function im1 = rotate_image(im2, R)
% im1 = rotate_image(im2, R)
% rotate image using rotation matrix
% im2 is the image to rotate
% R is the rotation matrix
% im1 is the rotated image

npixels = size(im2, 1);
dim = size(R, 1);

% compute corresponding rotation+translation from rotation matrix
if dim == 3
    [theta, dx, dy] = find_specs_3d(R);
else % dim == 2
    theta = find_specs_2d(R);
    dx = 0;
    dy = 0;
end

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);

% calculate corresponding affine transformation
A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; dx dy 1]);

% transform image
im1 = imwarp(im2, RI, A, 'outputview', RI);

% compute translation+rotation corresponding to a specific three-dimensional rotation matrix 
% R is the rotation matrix
% theta is the rotation angle
% dx is the x-translation (in fraction of the image)
% dy is the y-translation (in fraction of the image)
function [theta, dx, dy] = find_specs_3d(R)

angle_proj = 20;

alpha = atan2d(R(3,2), R(3,3));
gamma = atan2d(R(2,1), R(1,1));
beta = -asind(R(3,1));
% beta = atan2(-R(3,1), R(2,1)/sind(gamma));

theta = alpha;
dx = beta / angle_proj;
dy = gamma / angle_proj;

% compute rotation corresponding to a specific two-dimensional rotation matrix 
% R is the rotation matrix
% theta is the rotation angle
function theta = find_specs_2d(R)

theta = atan2d(R(2,1), R(2,2));



