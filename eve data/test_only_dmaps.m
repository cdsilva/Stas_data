clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';
print_figures = false;

eve_image_dir = '14_0216OreR_Eve_Dl';
%eve_membrane_times = 'Eve_13_0821/times_08_21_13.mat';

m = 60;

%%
% load(eve_membrane_times);
% t = times_08_21_13;

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);

image_channel = 2;

%figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', eve_image_dir, i));
    
    % resize image
    im1 = imresize(im1, [NaN npixels]);
    
    % extract relevent color from image
    im1 = im1(:,:,image_channel);
    
    npad = [npixels npixels]-size(im1);
    im1 = padarray(im1, floor(npad/2),'pre');
    im1 = padarray(im1, ceil(npad/2),'post');
        
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    %store image
    im1 = double(im1);
    image_set(:, :, i) = im1;
    
end

%% DMAPS

W = zeros(m);

for i=1:m
    for j=1:i-1
        W(i,j) = sum(sum((image_set(:,:,i) - image_set(:,:,j)).^2));
        W(j,i) = W(i,j);
    end
end

eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

%% order

[~, I] = sort(V(:,2));

figure;
for i=1:m 
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(uint8(image_set(:,:,I(i))));
    
end

figure;
for i=1:m-1
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(imabsdiff(uint8(image_set(:,:,I(i))),uint8(image_set(:,:,I(i+1)))));
    
end