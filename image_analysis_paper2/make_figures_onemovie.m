clear all
close all

%%

npixels = 100;

dt = 0.5;

nmovies = 6;

idx_end = [56 60 47 44 46 42];
idx_start = idx_end - 41;

nrot = 36;
nshifts = 0;
shift_max = 0.1;

neigs = 10;
alpha = 0;

channel = 1;

k_print_fig = 2;

theta_err = zeros(nmovies, 1);
rank_corr = zeros(nmovies, 1);

%%
for k=1:nmovies
    load(sprintf('movie%d.mat', k));
    idx = idx_start(k):idx_end(k);
    
    time = time * dt;
    
    nimages = length(time);
    
    image_set = zeros(npixels, npixels, nimages, 'uint8');
    
    for i=1:nimages
        image_set(:, :, i) = image_fn(images(:, :, channel, i), npixels);
    end
    
    %
    image_set = image_set(:, :, idx);
    time = time(idx) - time(idx(1));
    nimages = length(time);
    
    %
    rng(21*k);
    theta = 360*rand(nimages, 1);
    idx = randperm(nimages);
    
    for i=1:nimages
        image_set(:, :, i) = imrotate(image_set(:, :, i), theta(i), 'crop');
    end
    
    %
    image_set = image_set(:, :, idx);
    time = time(idx);
    theta = theta(idx);
    
    %
    [R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
    dim = size(R, 1) / nimages;
    
    %
    W2 = W.^2;
    eps = median(W2(:))/10;
    [R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);
    
    %
    theta_opt = zeros(size(theta));
    for i=1:nimages
        R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
        theta_opt(i) = atan2d(R_tmp(2,1), R_tmp(1,1));
    end
    fprintf('movie %d, average error = %2.2f \n', k, std(mod(theta - theta_opt, 360)))
    theta_err(k) =  std(mod(theta - theta_opt, 360));
    
    %
    if corr(time, embed_coord(:,1)) < 0
        embed_coord(:,1) = -embed_coord(:,1);
    end
    
    fprintf('movie %d, rank corr = %2.4f \n', k, corr(compute_ranks(time), compute_ranks(embed_coord(:,1))))
    rank_corr(k) = corr(compute_ranks(time), compute_ranks(embed_coord(:,1)));
    
    if k == k_print_fig
        [~, max_idx] = max(mod(theta, 360));
        make_fig(4.3,4);
        plot(mod(theta, 360), mod(theta_opt-theta_opt(max_idx)-1, 360)+1, '.')
        xlabel('true angle')
        ylabel('recovered angle')
        set(gca, 'xtick', [0 180 360]);
        set(gca, 'ytick', [0 180 360]);
        set(gca, 'xticklabel', {'0°'; '180°'; '360°'});
        set(gca, 'yticklabel', {'0°'; '180°'; '360°'});
        axis square
        saveas(gcf, 'angle_corr.pdf');
        
        make_fig(4,4);
        plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.')
        xlabel('true rank')
        ylabel('recovered rank')
        set(gca, 'xtick', [0 20 40])
        set(gca, 'ytick', [0 20 40])
        axis([0 50 0 50])
        axis square
        saveas(gcf, 'rank_corr.pdf');
        
        nprint_images = 10;
        subplot_dim1 = nprint_images;
        subplot_dim2 = 1;
        plot_idx = round(linspace(5, nimages, nprint_images));
        
        make_fig(17, 17/nprint_images);
        for j=1:nprint_images
            make_subplot(subplot_dim1, subplot_dim2, 0.01, j);
            imshow(image_set(:, :, plot_idx(j)));
        end
        saveas(gcf, 'movie_unregistered_unordered.pdf');
        
        make_fig(17, 17/nprint_images);
        [~, I] = sort(embed_coord(plot_idx, 1));
        for j=1:nprint_images
            im_tmp = imrotate(image_set(:, :, plot_idx(I(j))), -theta_opt(plot_idx(I(j)))+88, 'crop');
            make_subplot(subplot_dim1, subplot_dim2, 0.01, j);
            imshow(im_tmp);
            text(npixels/2, npixels/2, sprintf('%2.1f min', time(plot_idx(I(j)))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
        end
        saveas(gcf, 'movie_registered_ordered.pdf');
        return
    end
end

% make_fig(5.5, 5.5);
% bar(theta_err([2 1 3 5 6]));
% xlabel('movie index')
% ylabel('average angle error (degrees)')
% saveas(gcf, 'angle_corr_allmovies.pdf');
%
% make_fig(5.5, 5.5);
% bar(rank_corr([2 1 3 5 6]));
% xlabel('movie index')
% ylabel('rank correlation coefficient')
% saveas(gcf, 'rank_corr_allmovies.pdf');
%
% make_fig(5.5, 5.5)
% bar_width = 0.15;
% [ax,h1,h2] = plotyy((1:nmovies)-bar_width/2,theta_err,(1:nmovies)+bar_width/2,rank_corr,'bar');
% set(ax, 'xtick', 1:nmovies)
% set(h1, 'facecolor','r')
% set(h1, 'barwidth', 0.5)
% set(h2, 'barwidth', 0.5)

%%

load(sprintf('movie%d.mat', k_print_fig));
idx = idx_start(k_print_fig):idx_end(k_print_fig);

time = time * dt;

nimages = length(time);

image_set = zeros(npixels, npixels, nimages, 'uint8');

for i=1:nimages
    image_set(:, :, i) = image_fn(images(:, :, channel, i), npixels);
end

%
image_set = image_set(:, :, idx);
time = time(idx);
nimages = length(time);

%
rng(21*k);
theta = 360*rand(nimages, 1);
idx = randperm(nimages);

for i=1:nimages
    image_set(:, :, i) = imrotate(image_set(:, :, i), theta(i), 'crop');
end

%
image_set_all = image_set(:, :, idx);
time_all = time(idx);
nimages_all = length(time_all);
theta_all = theta(idx);

%
[R_all, W_all] = compute_pairwise_alignments(image_set_all, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages_all;

%
nbootstrap_samples = 200;
nbootstrap_points = 20;

sample_vec = round(linspace(10, nimages_all, nbootstrap_points));
theta_err = zeros(length(sample_vec), nbootstrap_samples);
rank_corr = zeros(length(sample_vec), nbootstrap_samples);

for j = 1:nbootstrap_points
    nimages = sample_vec(j);
    
    j
    
    for j2 = 1:nbootstrap_samples
        
        sample_idx = randi(nimages_all, 1,  nimages);
        
        image_set = image_set_all(:, :, sample_idx);
        time = time_all(sample_idx);
        theta = theta_all(sample_idx);
        
        %
        R_idx = reshape(repmat(dim*(sample_idx-1), dim, 1) + repmat((1:dim)', 1, nimages), [], 1);
        R = R_all(R_idx, R_idx);
        W = W_all(sample_idx, sample_idx);
        
        %
        W2 = W.^2;
        eps = median(W2(:))/10;
        [R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);
        
        %
        theta_opt = zeros(size(theta));
        for i=1:nimages
            R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
            theta_opt(i) = atan2d(R_tmp(2,1), R_tmp(1,1));
        end
        
        %
        theta_err(j, j2) = std(mod(theta-theta_opt, 360));
        if theta_err(j, j2) > 20
            theta_err(j, j2) = std(mod(theta-theta_opt+180, 360)-180);
        end
        
        %
        if corr(time, embed_coord(:,1), 'type','spearman') < 0
            embed_coord(:,1) = -embed_coord(:,1);
        end
        
        rank_corr(j, j2) = corr(compute_ranks(time), compute_ranks(embed_coord(:,1)));
    end
end

%
% figure;
% plot(sample_vec, mean(theta_err, 2), '.')
% xlabel('number of images')
% ylabel('average error in recovered angle')

make_fig(4,4);
plot(sample_vec, mean(rank_corr, 2), '.')
xlabel('number of images')
ylabel('rank correlation')
axis([0 50 0.5 1.0])
set(gca, 'xtick', [0 20 40])
axis square
saveas(gcf, 'bootstrap_rankcorr.pdf');

% figure;
% errorbar(sample_vec, mean(theta_err, 2), std(theta_err,[], 2))
%
% figure;
% errorbar(sample_vec, mean(rank_corr, 2), std(rank_corr,[], 2))
