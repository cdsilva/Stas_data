clear all
close all

%% directories where things are stored

im_save_dir = 'paper_figures2';

image_dir = '../membrane_pictures/14_0501_dpERK_late';

m = 132;
ind = 1:m;

%% image parameters

npixels = 100;

%set image plotting parameters
subplot_dim1 = floor(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim2:-1:1)/subplot_dim2, (1:subplot_dim1)/subplot_dim1);

nchannels = 3;
image_set = zeros(npixels, npixels, nchannels, m, 'uint8');
image_set_raw = zeros(1024, 1024, nchannels, m, 'uint8');

%% load images
figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    image_set_raw(:, :, :, i) = im1;
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    %     blur image
    H = fspecial('disk',5);
    im1(:, :, 1) = imfilter(im1(:, :, 1),H,'replicate');
    im1(:, :, 1) = 0.5*imadjust(im1(:, :, 1));
    
    im1(:, :, 3) = 1.5 * im1(:, :, 3);
    
    im1 = circshift(im1,[0 0 -1]);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    %subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    %store image
    image_set(:, :, :, i) = im1;
    
end



%% select which images to use

% ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 32 102 101 100]);
ind = 20:30;

image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);

m = length(ind);

% subplot_dim1 = floor(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);
subplot_dim1 = 18;
subplot_dim2 = 6;

[Y, X] = meshgrid((subplot_dim2:-1:1)/subplot_dim2, (1:subplot_dim1)/subplot_dim1);


%%

image_set_ordered = register_and_order_images(image_set);

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    
    imshow(image_set_ordered(:,:,:,i));
end

