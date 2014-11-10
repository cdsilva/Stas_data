% clear all
close all

% load fixed_images1.mat
images = squeeze(images2(:,:,1,:));


npixels = size(images, 1);

im1 = medfilt2(adapthisteq(images2(:,:,8)), [10 10]);
% im2 = medfilt2(adapthisteq(images_movie(:,:,10)), [10 10]);
im2 = medfilt2(images_movie(:,:,10), [10 10 ]);

[u, v] = HS(im1, im2);

[X, Y] = meshgrid(1:npixels,1:npixels);

%%
im2_warped = uint8(interp2(X, Y, double(im2), X+u, Y+v));

im_tmp = zeros(npixels, npixels, 3, 'uint8');
im_tmp(:,:,1) = im1;
im_tmp(:,:,2) = im2;

figure;
subplot(1,2,1)
imshow(im_tmp)

im_tmp(:,:,2) = im2_warped;
subplot(1,2,2)
imshow(im_tmp)

%%

navg = 20;
stride = npixels/navg;

avg_X = zeros(navg);
avg_Y = zeros(navg);
avg_u = zeros(navg);
avg_v = zeros(navg);
avg_im1 = zeros(navg);

for i=1:navg
    for j=1:navg
        avg_X(i,j) = mean(mean(X((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_Y(i,j) = mean(mean(Y((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_u(i,j) = mean(mean(u((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_v(i,j) = mean(mean(v((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_im1(i,j) = mean(mean(double(im1((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride))));
    end
end

%%

[~, I]= sort(avg_im1(:));
I = I(round(navg^2/2):end);
    
%%

sigma = 0.05;
% im1_warped = warp_image_rbf(im1, [reshape(avg_X, [], 1) reshape(avg_Y, [], 1)], [reshape(avg_X+avg_u, [], 1) reshape(avg_Y+avg_v, [], 1)], sigma);
im1_warped = warp_image_rbf(im1, [avg_X(I) avg_Y(I)], [avg_X(I)+avg_u(I) avg_Y(I)+avg_v(I)], sigma);

%%

im_tmp = zeros(npixels, npixels, 3, 'uint8');
im_tmp(:,:,1) = im1;
im_tmp(:,:,2) = im2;

figure;
subplot(1,2,1)
imshow(im_tmp)

im_tmp(:,:,1) = im1_warped;
subplot(1,2,2)
imshow(im_tmp)
