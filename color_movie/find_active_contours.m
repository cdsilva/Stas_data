function B = find_active_contours(im1)

npixels = 100;

im1 = imresize(im1, [npixels npixels]);
im1 = medfilt2(adapthisteq(im1), [3 3]);

mask = false(npixels);
[X, Y] = meshgrid(1:npixels, 1:npixels);
X = X - npixels/2;
Y = Y - npixels/2;
R2 = X.^2 + Y.^2;
idx = find(R2 < (npixels/2)^2 & R2 > (0.8*npixels/2)^2);
mask(idx) = true;

% mask = imadd(mask, im2bw(medfilt2(im1, [5 5]), 0.05));
% mask = im2bw(medfilt2(im1, [5 5]), 0.05);
smoothfactor = 2;

B = bwboundaries(activecontour(im1, mask, 'Chan-Vese', smoothfactor));