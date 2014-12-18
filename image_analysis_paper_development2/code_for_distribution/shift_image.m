function im2 = shift_image(im1)

npixels = size(im1, 1);
npixels2 = size(im1, 2);

bw = edge(rgb2gray(im1), 'canny');
bw = imfill(bw,'holes');
% bw = bwareaopen(bw,150);
bw = imclose(bw, strel('disk',round(0.01*(npixels+npixels2))));

[L,N] = bwlabel(bw,8);
imdata = regionprops(L,'Area','Centroid');
[~, i] = max([imdata.Area]);

tform = affine2d([1 0 0; 0 1 0; round(npixels/2)-imdata(i).Centroid(1) round(npixels2/2)-imdata(i).Centroid(2) 1]);

% im2 = imtranslate(im1, round(npixels/2)-imdata(i).Centroid(1), round(npixels2/2)-imdata(i).Centroid(2));
Rout = imref2d(size(im1));
im2 = imwarp(im1,tform, 'outputview', Rout, 'FillValues', squeeze(uint8(mean(mean(double(im1), 1), 2))));

