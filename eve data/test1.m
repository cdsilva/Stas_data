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
    
    %im1 = imrotate(im1, rand*360, 'crop');
    
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
neigs = 10;

%R_opt = ang_synch(R, dim);
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

toc

%%
image_set_aligned = zeros(size(image_set));
for i=1:m
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
    image_set_aligned(:,:,i) = double(im1);
end

%%
% W = zeros(m);
% for i=1:m
%     for j=1:i-1
%         W(i,j) = sum(sum((image_set_aligned(:,:,i)-image_set_aligned(:,:,j)).^2));
%         W(j,i) = W(i,j);
%     end
% end
% 
% eps = median(W(:));
% 
% [V, D] = dmaps(W, eps, 10);

%%

idx = find(embed_idx(1,:)==4 & embed_idx(2,:)==1);
%[~, idx ] = max(abs(corr(embed_coord, t)));

[~, I] = sort(embed_coord(:, idx));

% W = squareform(pdist(embed_coord)).^2;
% eps = median(W(:));
% [V, D] = dmaps(W, eps, 10);
% [~, I] = sort(V(:,3));

figure;
for i=1:m
    im1 = rotate_image(R_opt(dim*(I(i)-1)+1:dim*I(i),:), uint8(image_set(:,:,I(i))), angle_proj);
    %det(R_opt(dim*(i-1)+1:dim*i,:))
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end