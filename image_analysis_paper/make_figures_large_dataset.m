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
im_save_dir = 'paper_figures';

image_channel = 3;

%% load membrane lengths

load(dpERK_data);
mem_lengths = length;
clear length;

ind = setdiff(1:90, [3, 20, 25, 32, 53, 59, 61, 66, 82]);

m = length(ind);

mem_lengths = mem_lengths(ind);

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);


%figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir, ind(i)));
    
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

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 20;
shift_step = 4;
dim = 3;

tic
[R, W] = compute_pairwise_alignments(image_set, angle_proj, shift_max, shift_step);

%% VDM

eps = median(W(:));
neigs = 5;

%R_opt = ang_synch(R, dim);
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

toc

%%

figure;
for i=1:m
    im1 = uint8(image_set(:,:,i));
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end

image_set_aligned = zeros(size(image_set));
for i=1:m
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
    image_set_aligned(:,:,i) = double(im1);
end

%% order using vdm
[~, idx] = max(abs(corr(embed_coord, mem_lengths)));
%idx = 7;

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

[~, I] = sort(embed_coord(:, idx));
figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_set_aligned(:,:,I(i))));
end
if print_figures
    if image_channel == 2
        print(sprintf('%s/dpERK_array', im_save_dir), fmt, res);
    end
    if image_channel == 3
        print(sprintf('%s/DI_array', im_save_dir), fmt, res);
    end
end

if print_figures
    figure;
    set(gcf, 'paperposition',[0 0 8 8])
    for i=1:m
        imshow(uint8(image_set_aligned(:,:,I(i))));
        set(gca,'position',[0 0 1 1],'units','normalized')
        if image_channel == 2
            print(sprintf('%s/dpERK_%d', im_save_dir,i), fmt, res);
        end
        if image_channel == 3
            print(sprintf('%s/DI_%d', im_save_dir, i), fmt, res);
        end
        clf;
    end
end

r1 = compute_ranks(mem_lengths);
r2 = compute_ranks(embed_coord(:,idx));

disp(corr(r1, r2))

figure;
plot(r1, r2, '.')
xlabel('rank from membrane lengths')
if image_channel == 2
    ylabel('rank from dmaps using dpERK data')
    if print_figures
        print(sprintf('%s/dpERK_rank_corr', im_save_dir), fmt, res);
    end
end
if image_channel == 3
    ylabel('rank from dmaps using DI data')
    if print_figures
        print(sprintf('%s/DI_rank_corr', im_save_dir), fmt, res);
    end
end

%% order using dmaps on vdm

% W3 = squareform(pdist(embed_coord)).^2;
% eps3 = median(W3(:));
% 
% [V, D] = dmaps(W3, eps3, 10);
% 
% if corr(mem_lengths, V(:,2)) < 0
%     V(:,2) = -V(:,2);
% end
% 
% [~, I] = sort(V(:,2));
% figure;
% for i=1:m
%     subplot(subplot_dim1, subplot_dim2, i);
%     imshow(uint8(image_set_aligned(:,:,I(i))));
% end
% 
% r1 = compute_ranks(mem_lengths);
% r2 = compute_ranks(V(:,2));
% 
% disp(corr(r1, r2))

%% do PCA on VDM
% [V, D] = pca(embed_coord);
% 
% embed_coord2 = embed_coord - repmat(mean(embed_coord), m, 1);
% [~, I] = sort(embed_coord2*V(:,1));
% 
% r1 = compute_ranks(mem_lengths);
% r2 = compute_ranks(embed_coord2*V(:,1));
% 
% disp(corr(r1, r2))
% 
% figure;
% for i=1:m
%     subplot(subplot_dim1, subplot_dim2, i);
%     imshow(uint8(image_set_aligned(:,:,I(i))));
% end

%% order using dmaps

% W2 = zeros(m);
% for i=1:m
%     for j=1:i-1
%         W2(i,j) = sum(sum((image_set_aligned(:,:,i) - image_set_aligned(:,:,j)).^2));
%         W2(j,i) = W2(i,j);
%     end
% end
% eps2 = median(W2(:));
% [V, D] = dmaps(W2, eps2, 10);
%
% if corr(mem_lengths, V(:,2)) < 0
%    V(:,2) = -V(:,2);
% end
%
% [~, I] = sort(V(:,2));
% figure;
% for i=1:m
%     subplot(subplot_dim1, subplot_dim2, i);
%     imshow(uint8(image_set_aligned(:,:,I(i))));
% end
%
% r1 = compute_ranks(mem_lengths);
% r2 = compute_ranks(V(:,2));
%
% disp(corr(r1, r2))
%
% figure;
% plot(r1, r2, '.')
% xlabel('rank from membrane lengths')
% if image_channel == 2
%     ylabel('rank from dmaps using dpERK data')
% end
% if image_channel == 3
%     ylabel('rank from dmaps using DI data')
% end