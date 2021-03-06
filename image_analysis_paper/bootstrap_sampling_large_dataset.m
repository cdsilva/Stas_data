clear all
close all

res = '-r300';
fmt = '-djpeg';

im_save_dir = 'paper_figures';

%% read in data

% % directories where things are stored
% time_data = '../membrane_pictures/large_dataset/time.mat';
% image_dir = '../membrane_pictures/large_dataset';
% 
% % load membrane lengths
% load(time_data);
% mem_lengths = length;
% clear length;
% 
% % remove some younger embryos
% ind = setdiff(1:90, [3, 20, 25, 32, 53, 59, 61, 66, 82]);
% 
% m = length(ind);
% mem_lengths = mem_lengths(ind);
% 
% % compute ranks from membrane lengths
% ranks_from_membranes = compute_ranks(mem_lengths);
% 
% % image parameters
% npixels = 100;
% image_channel = 2;
% 
% %set image plotting parameters
% subplot_dim1 = ceil(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);
% 
% image_set = zeros(npixels, npixels, m);
% 
% %% load images
% for i=1:m
%     % read image
%     im1 = imread(sprintf('%s/lat%02d.tif', image_dir, ind(i)));
%     
%     % resize image
%     im1 = imresize(im1, [npixels npixels]);
%     
%     % extract relevent color from image
%     im1 = im1(:,:,image_channel);
%     
%     %store image
%     im1 = double(im1);
%     image_set(:, :, i) = im1;
%     
% end
% 
% % raw dpERK images-- synchronization
% angle_proj = pi/8;
% shift_max = 20;
% shift_step = 4;
% dim = 3;
% 
% [R, W] = compute_pairwise_alignments(image_set, angle_proj, shift_max, shift_step);

%% load
%load large_data_set.mat
load dpERK_aligned_all.mat

%% define number of bootstrap samples and bootstrap range

npoints_boot = 20;
boot_range = round(linspace(0.1*m, m, npoints_boot));

nsamples_boot = 50;

rng(12345);

r = zeros(npoints_boot, 1);

neigs = 5;

for i=1:npoints_boot
    for j=1:nsamples_boot
        
        % subsample data
        idx_boot = randsample(m, boot_range(i), false);
        idx_boot_R = reshape([dim*(idx_boot-1)+1 dim*(idx_boot-1)+2 idx_boot*dim]', [], 1);
        
        W2 = W(idx_boot, idx_boot);
        R2 = R(idx_boot_R, idx_boot_R);
        mem_lengths2 = mem_lengths(idx_boot);
        
        % VDM
        eps = median(W2(:));
        [R_opt, embed_coord, embed_idx, D] = vdm(R2, W2, eps, neigs);
        
        % order using vdm
        idx = find(embed_idx(1,:) == dim+1 & embed_idx(2,:) == 1);
        
        if corr(mem_lengths2, embed_coord(:, idx)) < 0
            embed_coord(:, idx) = -embed_coord(:, idx);
        end
        
        r1 = compute_ranks(mem_lengths2);
        r2 = compute_ranks(embed_coord(:,idx));
        
        r(i) = r(i) + corr(r1, r2);
    end
end

r = r / nsamples_boot;

figure;
plot(boot_range, r, '.')
xlabel('number of data points')
ylabel('average rank correlation')
% if print_figures
%     print(sprintf('%s/bootstrap_rank_corr', im_save_dir), fmt, res);
% end

