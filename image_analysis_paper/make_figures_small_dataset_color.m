clear all
close all

res = '-r300';
fmt = '-djpeg';

fontsize = 10;

im_save_dir = 'paper_figures';

%% directories where things are stored
time_data = '../membrane_lengths/oct16.mat';
image_dir = '../membrane_pictures/membrane2/dpERK_staining';

%% load membrane lengths

load(time_data);
mem_lengths = L(:,1);

% remove some embryos
ind = setdiff(1:52, [9, 15, 22]);

m = length(ind);

mem_lengths = mem_lengths(ind);

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% image parameters

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim1:-1:1)/subplot_dim1, (1:subplot_dim2)/subplot_dim2);

nchannels = 3;
image_set = zeros(npixels, npixels, nchannels, m, 'uint8');

nuclei = zeros(npixels, npixels, m, 'uint8');

%% load images

for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    if i==1
        im1 = circshift(im1, [0 15 0]);
    end
    
    nuclei(:, :, i) = im1(:, :, 1);
    
    im1(:, :, 1) = im1(:, :, 2);
    im1(:, :, 2) = im1(:, :, 3);
    im1(:, :, 3) = 0;
    
    %store image
    image_set(:, :, :, i) = im1;
    
end

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 20;
shift_step = 4;
dim = 3;

% [R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% save('pairwise_alignments_small.mat', 'R', 'W');
load('pairwise_alignments_small.mat');

%% synchronization

R_opt = ang_synch(R, dim);

%%

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set(:,:,:,i);
    
    for j=1:3
        im1(:,:,j) = im1(:,:,j) + nuclei(:,:,i);
    end
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end
%print(sprintf('%s/unregistered_unordered_2d', im_save_dir), fmt, res);
saveas(gcf, sprintf('%s/fig2a', im_save_dir), 'epsc');
saveas(gcf, sprintf('%s/fig2a', im_save_dir), 'tif');

image_set_aligned = zeros(size(image_set), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    for j=1:3
        im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set_aligned(:,:,:,i);
    
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end
%print(sprintf('%s/registered_unordered_2d', im_save_dir), fmt, res);

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

eps2 = median(W2(:));
neigs = 5;

[V, D] = dmaps(W, eps2, neigs);

if corr(mem_lengths, V(:,2)) < 0
    V(:,2) = -V(:,2);
end

[~, I] = sort(V(:,2));
figure;
set(gcf, 'paperposition', [0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1 1/subplot_dim2])
    
    imshow(image_set_aligned(:,:,:,I(i)));
end
%print(sprintf('%s/registered_ordered_2d', im_save_dir), fmt, res);


ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_dmaps = compute_ranks(V(:,2));

figure;
plot(ranks_from_membranes, ranks_from_dmaps, '.')
xlabel('rank from membrane lengths')
ylabel('rank from synchronization + dmaps')
%print(sprintf('%s/rank_corr_dmaps', im_save_dir), fmt, res);


%% VDM

eps = median(W(:));

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

image_set_aligned = zeros(size(image_set), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    for j=1:3
        im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
end

%% order using vdm

idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

[~, I] = sort(embed_coord(:, idx));
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    
    imshow(image_set_aligned(:,:,:,I(i)));
end
%print(sprintf('%s/registered_ordered_small_vdm_2d', im_save_dir), fmt, res);
saveas(gcf, sprintf('%s/fig5a', im_save_dir), 'epsc');
saveas(gcf, sprintf('%s/fig5a', im_save_dir), 'tif');

ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_vdm = compute_ranks(embed_coord(:,idx));

figure;
set(gcf, 'papersize', [11/4 8.5/4])
set(gcf, 'paperposition', [0 0 11/4 8.5/4])
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('rank from cellularization')
ylabel('rank from vdm')
%print(sprintf('%s/rank_corr_vdm', im_save_dir), fmt, res);

%%
% 
% for i=1:m
%     im_tmp = imread(sprintf('%s/emb%02d.tif', image_dir, ind(I(i))));
% 
%     im1 = rotate_image(R_opt(dim*(I(i)-1)+1:dim*I(i),:), im_tmp, angle_proj);
% 
%     imwrite(im1, sprintf('%s/registered_emb_%02d.tif', im_save_dir, i));
% end

%%

load('../membrane_pictures/large_dataset/dpERK_DL_140313.mat');
dpERK = dpERK_oct16;
Dl = Dl_oct16;

dpERK = dpERK(ind,:);
Dl = Dl(ind, :);

[m, n] = size(dpERK);

data_1d = zeros(m, n, 3);
data_1d(:,:,1) = dpERK;
data_1d(:,:,2) = Dl;

for i=1:2
    data_1d(:,:,i) = data_1d(:,:,i) / max(max(data_1d(:,:,i))) * 256;
end

figure;
imshow(uint8(data_1d))

[~, I] = sort(mem_lengths);
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
image(uint8(data_1d(I, :, :)))
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('position', 'fontsize', fontsize)
ylabel('sample', 'fontsize', fontsize)
saveas(gcf, sprintf('%s/fig6b', im_save_dir), 'epsc');

%%
%test_idx = 41;
test_idx = 27;

im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(test_idx)));
nuclei_test = im1(:, :, 1);

im1(:, :, 1) = im1(:, :, 2);
im1(:, :, 2) = im1(:, :, 3);
im1(:, :, 3) = 0;

for i=1:3
    im1(:, :, i) = im1(:, :, i) + nuclei_test;
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
imshow(im1);
saveas(gcf, sprintf('%s/fig2b', im_save_dir), 'epsc');

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
plot(data_1d(test_idx, :, 1), '-r')
hold on
plot(data_1d(test_idx, :, 2), 'color', [0.2 0.9 0.2])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
ylabel('intensity', 'fontsize', fontsize)
xlabel('position', 'fontsize', fontsize)
saveas(gcf, sprintf('%s/fig2c', im_save_dir), 'epsc');



