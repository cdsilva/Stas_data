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

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);

image_channel = 2;

figure;
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
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    image_set2(:,:,i) = double(im1);
    
end

%% raw images-- compute scattering transform coefficients

addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
%addpath '../../../MATLAB/scatnet-0.2';
addpath_scatnet

% set scattering transform parameters
filt_opt = struct();
filt_rot_opt = struct();
% oversampling is set to infty
% so that scattering of original and rotated
% image will be sampled at exactly the same points
scat_opt = struct();
scat_opt.oversampling = 10;

% compute scattering transform of first image
x = image_set(:,:,1);
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
sx_all = zeros(m, size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:m
    i
    x = image_set2(:,:,i);
    
    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%% dmaps
W = squareform(pdist(sx_all)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

% if corr(V(:,2), L(:,1)) < 0
%     V(:,2) = -V(:,2);
% end

[~, idx] = max(abs(corr(V, t)));
[~, I] = sort(V(:,idx));

figure;
for i=1:m
    im1 = uint8(image_set2(:,:,I(i)));
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end
