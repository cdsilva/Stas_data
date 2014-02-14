clear all
close all

dpERK_data = '../membrane_pictures/large_dataset/time.mat';
dpERK_image_dir = '../membrane_pictures/large_dataset';

%% load membrane lengths

load(dpERK_data);
mem_lengths = length;
clear length;

m   = 16;

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);

image_channel = 2;

%%
[~, ind] = max(mem_lengths);

im_template = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir,ind));
%resize image
im_template = imresize(im_template, [npixels npixels]);
%extract relevent color from image
im_template = im_template(:,:,image_channel);

theta_vec = linspace(0, 360, 21);
theta_vec = theta_vec(1:end-1);

thetas = zeros(m,1);
xdisp = zeros(m,1);
ydisp = zeros(m,1);

figure;
for i=1:m
    thetas(i) = theta_vec(randi(20,1));
    xdisp(i) = randi([-10 10], 1);
    ydisp(i) = randi([-10 10], 1);
    
    im1 = imrotate(im_template, thetas(i), 'crop');
    im1 = circshift(im1, [xdisp(i) ydisp(i)]);
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    %store image
    im1 = double(im1);
    image_set(:, :, i) = im1;
    
end

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 20;
shift_step = 4;
dim = 3;

tic
[R, W] = compute_pairwise_alignments(image_set, angle_proj, shift_max, shift_step);

%%

eps = median(W(:));
neigs = 5;

%R_opt = ang_synch(R, dim);
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

toc


figure;
for i=1:m
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
    %det(R_opt(dim*(i-1)+1:dim*i,:))
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end