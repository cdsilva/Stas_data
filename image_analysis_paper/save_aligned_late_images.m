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
    
    %store image
    image_set(:, :, :, i) = im1;
    
end

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 10;
shift_step = 2;
dim = 3;

% matlabpool open 4;
% [R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% matlabpool close;
%
%save('pairwise_alignments_late_noblur.mat', 'R', 'W');
%save('pairwise_alignments_late.mat', 'R', 'W');
%save('pairwise_alignments_nonuclei.mat', 'R', 'W');



load('pairwise_alignments_late.mat');

%load('pairwise_alignments_nonuclei.mat');


%% select which images to use

ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 32 102 101 100]);

image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);

W = W(ind, ind);

R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
R = R(R_ind, R_ind);

m = length(ind);

% subplot_dim1 = floor(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);
subplot_dim1 = 18;
subplot_dim2 = 6;

[Y, X] = meshgrid((subplot_dim2:-1:1)/subplot_dim2, (1:subplot_dim1)/subplot_dim1);

%% VDM

eps = median(W(:))/10;
neigs = 3*m;
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

image_set_aligned = zeros(size(image_set), 'uint8');
image_set_raw_aligned = zeros(size(image_set_raw), 'uint8');
theta_adjust = 90;

for i=1:m
    im_tmp = image_set(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    im1 = circshift(im1, [5 10 0]);
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set_raw(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    im1 = circshift(im1, [50 100 0]);
    
    image_set_raw_aligned(:,:,:,i) = im1;
end

%% order using vdm

idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);

[~, I] = sort(embed_coord(:, idx));
if find(I == 10) < m/2
    I = flipud(I);
end

image_set_raw_aligned_ordered = image_set_raw_aligned(:,:,:,I);

save('late_images_aligned.mat','image_set_raw_aligned_ordered');