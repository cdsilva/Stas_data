clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';
print_figures = false;

% eve_image_dir = 'Eve_13_0821';
% eve_membrane_times = 'Eve_13_0821/times_08_21_13.mat';
% load(eve_membrane_times);
% t = times_08_21_13;

eve_image_dir = '14_0216OreR_Eve_Dl';
eve_membrane_times = '14_0216OreR_Eve_Dl/times_14_02_16.mat';
load(eve_membrane_times);
t = times_14_02_16;

m = length(t);

%% load images

npixels = 60;

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

%% PCA

pca_data = zeros(m, npixels^2);

for i=1:m
    pca_data(i,:) = reshape(image_set(:,:,i), 1, []);
end

pca_data = pca_data - repmat(mean(pca_data), m, 1);



[V, D] = PCA(pca_data, 5);

figure;
for i=1:4
    subplot(2, 2, i)
    imagesc(reshape(V(:,i), npixels, []));
    colormap(gray)
    axis off
end

figure;
for i=1:4
    subplot(2,2, i)
    imagesc(abs(reshape(V(:,i), npixels, [])));
    colormap(gray)
    axis off
end

%%
[X, Y] = meshgrid(1:npixels,1:npixels);
X = (X-mean(1:npixels))/(npixels/2);
Y = (Y-mean(1:npixels))/(npixels/2);

idx = find(X.^2/0.5+Y.^2/0.1 > 1);

image_set2 = zeros(size(image_set));

figure;
for i=1:m
    im1 = uint8(image_set(:, :, i));
    
    im1(idx) = 0;
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    image_set2(:,:,i) = double(im1);
    
end

%%

pca_data = zeros(m, npixels^2);

for i=1:m
    pca_data(i,:) = reshape(image_set2(:,:,i), 1, []);
end
pca_data = pca_data - repmat(mean(pca_data), m, 1);


[V, D] = PCA(pca_data, 5);

figure;
for i=1:4
    subplot(2, 2, i)
    imagesc(reshape(V(:,i), npixels, []));
    colormap(gray)
    axis off
end

figure;
for i=1:4
    subplot(2,2, i)
    imagesc(abs(reshape(V(:,i), npixels, [])));
    colormap(gray)
    axis off
end

figure; 
plot(t, pca_data*V(:,1),'.')



