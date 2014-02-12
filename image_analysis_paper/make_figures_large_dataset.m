clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';
print_figures = false;

dpERK_data = '../membrane_pictures/large_dataset/time.mat';
dpERK_image_dir = '../membrane_pictures/large_dataset';
%dpERK_membrane_dir = '../membrane_pictures/membrane2';

%% load membrane lengths

load(dpERK_data);
mem_lengths = length;
clear length;

m = length(mem_lengths);
m = 9;

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);

image_channel = 2;

% image indices to plot
% im_save_idx = [1,19,25,30,35,40,45,49];

% figure;
% for i=1:m
%     % read image
%     im1 = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir, i));
% 
%     % resize image
%     im1 = imresize(im1, [npixels npixels]);
% 
%     % extract relevent color from image
%     im1 = im1(:,:,image_channel);
% 
%     subplot(subplot_dim1, subplot_dim2, i)
%     imshow(im1);
%     
%     %store image
%     im1 = double(im1);
%     image_set(:, :, i) = im1;
%     
% end


[~, ind] = max(mem_lengths);

im_template = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir,ind));
% resize image
im_template = imresize(im_template, [npixels npixels]);
% extract relevent color from image
im_template = im_template(:,:,image_channel);

theta_vec = linspace(0, 360, 21);
theta_vec = theta_vec(1:end-1);
figure;
for i=1:m
    
    im1 = imrotate(im_template, theta_vec(randi(20,1)), 'crop');
    im1 = circshift(im1, randi([-10 10], [1 2]));
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    %store image
    im1 = double(im1);
    image_set(:, :, i) = im1;
    
end

%%
% [~, ind] = sort(mem_lengths);
% ind = ind(end-10:end);
% 
% mem_lengths = mem_lengths(ind);
% image_set = image_set(:,:,ind);
% m = length(mem_lengths);
% 
% figure;
% for i=1:m
%     subplot(subplot_dim1, subplot_dim2, i)
%     imshow(uint8(image_set(:,:,i)));
% end

L = mem_lengths;

%% raw dpERK images-- synchronization

xmax = 0.1;
ymax = 0.1;
shift_max = 10;
shift_step = 2;
dim = 3;
[R, W] = compute_pairwise_alignments(image_set, xmax, ymax, shift_max, shift_step);

eps = median(W(:));
neigs = 5;
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

figure;
for i=1:m
    im1 = rotate_image(R_opt(1:dim,1:dim)'*R_opt(dim*(i-1)+1:dim*i,:), image_set(:,:,i), xmax, ymax);
    det(R_opt(1:dim,1:dim)'*R_opt(dim*(i-1)+1:dim*i,:))
    %im1 = rotate_image(R(dim*(i-1)+1:dim*i,4:6), image_set(:,:,i), xmax, ymax);
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(im1))
end

return
%make_2d_alignment_figures;

