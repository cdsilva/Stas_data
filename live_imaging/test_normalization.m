clear all 
close all

im1_raw = imread('../membrane_pictures/14_0501_dpERK_late/emb06.tif');

npixels_movie = 500;
[images, times] = read_video('bomyi_emb01_gast01.avi', npixels_movie);

%%
im2_raw = images(:,:,:,42);

channel = 1;
im1 = im1_raw(:,:,channel);
im2 = im2_raw(:,:,channel);

figure;
subplot(1,2,1)
imshow(im1);
subplot(1,2,2)
imshow(im2);


% H = fspecial('disk', 3);
% image_fn = @(image) adapthisteq(imfilter(imresize(image, [100 100]), H, 'replicate'));

% H = fspecial('disk', 25);
% image_fn = @(image) adapthisteq(imfilter(imresize(image, [1000 1000]), H, 'replicate'),'numtiles', [4 4], 'cliplimit', 0.01);

H = fspecial('disk', 25);
image_fn = @(image) adapthisteq(imfilter(imresize(image, [1000 1000]), H, 'replicate'),'distribution','exponential');

figure;
subplot(1,2,1)
imshow(image_fn(im1));
subplot(1,2,2)
imshow(image_fn(im2));
