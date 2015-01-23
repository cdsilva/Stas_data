function im1 = rotate_image(im2, R)
% im1 = rotate_image(im2, R)
% rotate image using rotation matrix
% im2 is the image to rotate
% R is the rotation matrix
% im1 is the rotated image

npixels = size(im2, 1);

theta = find_specs_2d(R);
dx = 0;
dy = 0;

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);

% calculate corresponding affine transformation
A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; dx dy 1]);

% transform image
im1 = imwarp(im2, RI, A, 'outputview', RI);

% compute rotation corresponding to a specific two-dimensional rotation matrix
% R is the rotation matrix
% theta is the rotation angle
function theta = find_specs_2d(R)

theta = atan2d(R(2,1), R(2,2));



