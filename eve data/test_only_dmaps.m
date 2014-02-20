clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';

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

[X, Y] = meshgrid(1:npixels,1:npixels);
X = (X-mean(1:npixels))/(npixels/2);
Y = (Y-mean(1:npixels))/(npixels/2);

idx = find(X.^2/0.5+Y.^2/0.1 > 1);

image_set2 = zeros(size(image_set));

figure;
for i=1:m
    im1 = uint8(image_set(:, :, i));
    
    im1(idx) = 0;
    
    im1 = imadjust(im1);
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    image_set2(:,:,i) = double(im1);
    
end


%% DMAPS

W = zeros(m);

for i=1:m
    for j=1:i-1
        W(i,j) = sum(sum((image_set2(:,:,i) - image_set2(:,:,j)).^2));
        W(j,i) = W(i,j);
    end
end

eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

%%

if corr(V(:,2), t) < 0
    V(:,2) = -V(:,2);
end

[~, I] = sort(V(:,2));

figure;
for i=1:m
    im1 = uint8(image_set2(:,:,I(i)));
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end
print('adjusted_images', fmt, res);

figure;
for i=1:m
    im1 = uint8(image_set(:,:,I(i)));
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end
print('raw_images', fmt, res);

figure;
plot(t, V(:,2),'.')

r1 = compute_ranks(t);
r2 = compute_ranks(V(:,2));

figure;
plot(r1, r2, '.')
xlabel('rank from membrane length')
ylabel('rank from dmaps')
print('rank_corr', fmt, res);

disp(corr(r1, r2))
