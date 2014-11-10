function im2 = warp_image_rbf2(im1, cp1, cp_deltas, sigma)
% this function warps an image using selected control points
% im1 is the final image
% cp1 is the location of the control points in the original image 
% (number of control points x 2 array)
% cp2 is the location of the control points in the warped image
% (number of control points x 2 array)
% sigma is the standard devation of the kernel used for the interpolation
% (sigma is in fraction of the image, e.g., sigma=0.1 means that the kernel
% spans about 10% of the image)
% im2 is the back-warped image

%% initialize parameters
% number of pixels in one dimension (assume square images)
npixels = size(im1, 1);
% number of control points
ncp = size(cp1, 1);

% locations of pixels in image
[X, Y] = meshgrid(1:npixels, 1:npixels);

%% compute rbf weights (use Gaussian kernel) for each control point at each pixel
weights = zeros(npixels, npixels, ncp);
for i=1:ncp
    % distance between control points and pixels
    weights(:,:,i) = (X - cp1(i,1)).^2 + (Y - cp1(i,2)).^2;
end
% Gaussian kernel
weights = exp(-weights/(npixels^2*sigma^2));

% normalize weights (weights for each pixel should sum to 1)
sum_weights = sum(weights, 3);
for i=1:ncp
    weights(:,:,i) = weights(:,:,i) ./ sum_weights;
end

%% compute shift at each pixel using control point shifts + rbf weights
X_deltas = zeros(size(X));
Y_deltas = zeros(size(Y));
for i=1:ncp
    X_deltas = X_deltas + weights(:,:,i) * cp_deltas(i,1);
    Y_deltas = Y_deltas + weights(:,:,i) * cp_deltas(i,2);
end

%% shift pixels, and interpolate back to original grid to get new image
if ndims(im1) == 2
    im2 = uint8(interp2(X, Y, double(im1), X+X_deltas,Y+Y_deltas));
else
    im2 = zeros(size(im1), 'uint8');
    for i=1:size(im1, 3)
        im2(:,:,i) = uint8(interp2(X,Y,double(im1(:,:,i)), X+X_deltas,Y+Y_deltas));
    end
end




