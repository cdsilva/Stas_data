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
subplot_dim1 = ceil(sqrt(m));
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
    %used on 5/5
    H = fspecial('gaussian',[5 5], 2.0);
    im1(:, :, 1) = imfilter(im1(:, :, 1),H,'replicate');
    im1(:, :, 1) = 0.5*imadjust(im1(:, :, 1));
    %     im1(:, :, 1) = 0;
    
    im1(:, :, 3) = 1.5 * im1(:, :, 3);
    
    im1 = circshift(im1,[0 0 -1]);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    %subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
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
load('pairwise_alignments_late.mat');


%% select which images to use
% all
%ind = 1:m;

% used for pictures on 5/6
ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 50 102]);

% used for pictures on 5/5
%ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105]);

% only young embryos
%ind = [55,29,123,35,85,19,54,58,93,69,74,17,105,7,70,14,131,115,102,67,82,112,97,117,81,40,127,63,124,96,104,48,120,121,5,92,111,41,34,1,3,2,9,130,20,21,83,108,76,60,44,62,71,30,26,49,36,64,12,57];

image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);

W = W(ind, ind);

R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
R = R(R_ind, R_ind);

m = length(ind);

subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim2:-1:1)/subplot_dim2, (1:subplot_dim1)/subplot_dim1);

%% synchronization

R_opt = ang_synch(R, dim);

%%

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    im1 = image_set(:,:,:,i);
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
end
saveas(gcf,sprintf('%s/raw_data2', im_save_dir), 'pdf')


image_set_aligned = zeros(size(image_set), 'uint8');
image_set_raw_aligned = zeros(size(image_set_raw), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set_raw(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_raw_aligned(:,:,:,i) = im1;
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set_raw_aligned(:,:,:,i);
    
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end

%% order

W2 = zeros(m);
for i=1:m
    imagei = image_set_aligned(:,:,:,i);
    for j=1:i-1
        imagej = image_set_aligned(:,:,:,j);
        d2 = sum((double(imagei(:))-double(imagej(:))).^2);
        W2(i, j) = d2;
        W2(j, i) = W2(i, j);
    end
end

eps2 = median(W2(:))/10;
neigs = 5;

[V, D] = dmaps(W, eps2, neigs);

[~, I] = sort(V(:,2));
figure;
set(gcf, 'paperposition', [0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1 1/subplot_dim2])
    
    imshow(image_set_raw_aligned(:,:,:,I(i)));
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1 1/subplot_dim2])
    
    imshow(image_set_aligned(:,:,:,I(i)));
end

%% VDM

eps = median(W(:))/10;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

image_set_aligned = zeros(size(image_set), 'uint8');
image_set_raw_aligned = zeros(size(image_set_raw), 'uint8');

for i=1:m
    im_tmp = image_set(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set_raw(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_raw_aligned(:,:,:,i) = im1;
end

%% order using vdm

idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);

[~, I] = sort(embed_coord(:, idx));
if find(I == 10) < m/2
    I = flipud(I);
end
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    
    imshow(image_set_raw_aligned(:,:,:,I(i)));
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.005 1/subplot_dim2-0.005])
    im1 = image_set_aligned(:,:,:,I(i));
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    
    imshow(im1);
end
saveas(gcf,sprintf('%s/VDM_ordered', im_save_dir), 'pdf')

%%

nstages = 12;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
    stage_indices = I(max(1, round((i-1)*m/nstages)):min(m,round(i*m/nstages)));
    im1 = uint8(mean(double(image_set_aligned(:,:,:,stage_indices)), 4));
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
    
end
saveas(gcf,sprintf('%s/average_trajectory', im_save_dir), 'pdf')


