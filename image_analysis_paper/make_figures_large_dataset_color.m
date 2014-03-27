clear all
close all

res = '-r300';
fmt = '-djpeg';

im_save_dir = 'paper_figures';

%% directories where things are stored
time_data = '../membrane_pictures/large_dataset/time.mat';
image_dir = '../membrane_pictures/large_dataset';

%% load membrane lengths

load(time_data);
mem_lengths = length;
clear length;

% remove some younger embryos
ind = setdiff(1:90, [3, 20, 25, 32, 53, 59, 61, 66, 82]);

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


%% load images

for i=1:m
    % read image
    im1 = imread(sprintf('%s/lat%02d.tif', image_dir, ind(i)));
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
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

%[R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
%save('pairwise_alignments.mat', 'R', 'W');
load('pairwise_alignments.mat');

%% synchronization

R_opt = ang_synch(R, dim);

%%

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set(:,:,:,i);
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end
print(sprintf('%s/unregistered_unordered_2d', im_save_dir), fmt, res);

image_set_aligned = zeros(size(image_set), 'uint8');
for i=1:m
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), image_set(:,:,:,i), angle_proj);
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
print(sprintf('%s/registered_unordered_2d', im_save_dir), fmt, res);

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
print(sprintf('%s/registered_ordered_2d', im_save_dir), fmt, res);


ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_dmaps = compute_ranks(V(:,2));

figure;
plot(ranks_from_membranes, ranks_from_dmaps, '.')
xlabel('rank from membrane lengths')
ylabel('rank from synchronization + dmaps')
print(sprintf('%s/rank_corr_dmaps', im_save_dir), fmt, res);

%save('all_aligned.mat')


%% VDM

eps = median(W(:));


[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
%% order using vdm

idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

[~, I] = sort(embed_coord(:, idx));
figure;
set(gcf, 'paperposition', [0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(image_set_aligned(:,:,:,I(i)));
end
print(sprintf('%s/registered_ordered_vdm_2d', im_save_dir), fmt, res);

ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_vdm = compute_ranks(embed_coord(:,idx));

figure;
set(gcf, 'papersize', [11/4 8.5/4])
set(gcf, 'paperposition', [0 0 11/4 8.5/4])
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('rank from cellularization')
ylabel('rank from vdm')
print(sprintf('%s/rank_corr_vdm', im_save_dir), fmt, res);

%save('all_aligned.mat')
