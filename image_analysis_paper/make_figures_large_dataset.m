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
%m = 16;

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);

image_channel = 3;

% image indices to plot
% im_save_idx = [1,19,25,30,35,40,45,49];

%figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir, i));

    % resize image
    im1 = imresize(im1, [npixels npixels]);

    % extract relevent color from image
    im1 = im1(:,:,image_channel);

    %subplot(subplot_dim1, subplot_dim2, i)
    %imshow(im1);
    
    %store image
    im1 = double(im1);
    image_set(:, :, i) = im1;
    
end

%%
% [~, ind] = max(mem_lengths);
% 
% im_template = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir,ind));
% %resize image
% im_template = imresize(im_template, [npixels npixels]);
% %extract relevent color from image
% im_template = im_template(:,:,image_channel);
% 
% theta_vec = linspace(0, 360, 21);
% theta_vec = theta_vec(1:end-1);
% 
% thetas = zeros(m,1);
% xdisp = zeros(m,1);
% ydisp = zeros(m,1);
% 
% figure;
% for i=1:m
%     thetas(i) = theta_vec(randi(20,1));
%     xdisp(i) = randi([-10 10], 1);
%     ydisp(i) = randi([-10 10], 1);
%     
%     im1 = imrotate(im_template, thetas(i), 'crop');
%     im1 = circshift(im1, [xdisp(i) ydisp(i)]);
%     
%     subplot(subplot_dim1, subplot_dim2, i)
%     imshow(im1);
%     
%     %store image
%     im1 = double(im1);
%     image_set(:, :, i) = im1;
%     
% end
% 

%%
m = 81;
[~, ind] = sort(mem_lengths);
ind = ind(end-m+1:end);
ind = ind(randperm(m));

mem_lengths = mem_lengths(ind);
image_set = image_set(:,:,ind);

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

% figure;
% for i=1:m
%     subplot(subplot_dim1, subplot_dim2, i)
%     imshow(uint8(image_set(:,:,i)));
% end

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

%%


[~, idx] = max(abs(corr(embed_coord, mem_lengths)));
idx = 7;

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

figure; 
plot(mem_lengths, embed_coord(:, idx),'.')
xlabel('membrane length')
ylabel(sprintf('$$ \\langle \\phi_%d, \\phi_%d \\rangle $$ (using DI data)', embed_idx(1, idx), embed_idx(2, idx)),'interpreter','latex')

[~, I] = sort(embed_coord(:, idx));
figure;
for i=1:m
    im1 = rotate_image(R_opt(dim*(I(i)-1)+1:dim*I(i),:), uint8(image_set(:,:,I(i))), angle_proj);
    %det(R_opt(dim*(i-1)+1:dim*i,:))
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end

figure;
for i=1:m
    im1 = uint8(image_set(:,:,i));
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end

%%
r1 = compute_ranks(mem_lengths);
r2 = compute_ranks(embed_coord(:,idx));

figure; 
plot(r1, r2, '.')
xlabel('rank from membrane lengths')
ylabel('rank from dmaps using DI data')