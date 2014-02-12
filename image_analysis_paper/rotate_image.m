function im1 = rotate_image(R, im2, xmax, ymax)

npixels = size(im2, 1);

[theta, dxpixels, dypixels] = find_specs(R, xmax, ymax, npixels);

im1 = circshift(im2, [dxpixels dypixels]);
im1 = imrotate(im1, theta, 'crop');


function [theta, dxpixels, dypixels] = find_specs(R, xmax, ymax, npixels)

% find rotation of delta x delta square projected on unit sphere
delta = 0.1;

% radius of sphere
r = 1;

% construct square
points = [-delta -delta;
    -delta delta;
    delta delta;
    delta -delta];

% project points onto unit sphere
points = [points sqrt(r^2-points(:,1).^2-points(:,2).^2)];

points_old = points;

% rotate points
points = (R*points')';

points = points(:,1:2);
points_old = points_old(:,1:2);

disp = mean(points) - mean(points_old);

dxpixels = round(disp(1)/(2*xmax) * npixels);
dypixels = round(disp(2)/(2*ymax) * npixels);

points = points - repmat(mean(points), size(points,1),1);

% calculate covariance matrix
C = points_old'*points;

% do svd
[u, s, v] = svd(C);

d = sign(det(C));
I = eye(2);
I(end,end) = d;

R = v * I * u';

theta = atan2(R(2,1),R(1,1))*(180/pi);


