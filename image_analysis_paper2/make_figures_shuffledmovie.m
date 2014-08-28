clear all
close all

%%

npixels = 100;

dt = 0.5;

movies_to_use = [2 1 3 5 6];
nmovies = 5;

nimages = 40;
idx_end = [56 60 47 44 46 42];
idx_start = idx_end - nimages + 1;

image_set = zeros(npixels, npixels, nimages, nmovies, 'uint8');
time_set = zeros(nimages, nmovies);

nrot = 36;
nshifts = 0;
shift_max = 0.1;

neigs = 10;
alpha = 0;

channel = 1;

%%

for k1 = 1:length(movies_to_use)
    k = movies_to_use(k1);
    
    load(sprintf('movie%d.mat', k));
    idx = idx_start(k):idx_end(k);
    
    images = images(:, :, :, idx);
    
    for i=1:nimages
        image_set(:, :, i, k1) = imrotate(image_fn(images(:, :, channel, i), npixels), 360*rand, 'crop');
    end
    
    time_set(:,k1) = dt * (time(idx) - time(idx(1)));
    
end

%%
figure;
for i=1:nmovies
    for j=1:nimages/4
        subplot(nmovies, nimages/4, (i-1)*nimages/4+j)
        imshow(image_set(:, :, 4*j, i))
    end
end


%%

nsubsamples = 200;
rank_corr = zeros(nsubsamples, 1);

for j=1:nsubsamples
    
    j
    
    % rng(333);
    rng(j);
    idx = randi(nmovies, nimages, 1);
    image_set_subsampled = zeros(npixels, npixels, nimages, 'uint8');
    time_subsampled = time_set(:,1);
    for i=1:nimages
        image_set_subsampled(:, :, i) = image_set(:, :, i, idx(i));
    end
    
    % figure;
    % for i=1:nimages
    %     subplot(6, 7,i)
    %     imshow(image_set_subsampled(:,:,i))
    % end
    
    %
    
    display_waitbar = false;
    [R, W] = compute_pairwise_alignments(image_set_subsampled, nrot, nshifts, shift_max, display_waitbar);
    dim = size(R, 1) / nimages;
    
    %
    W2 = W.^2;
    eps = median(W2(:))/10;
    [R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);
    
    %
    image_set_aligned = zeros(npixels, npixels, nimages, 'uint8');
    for i=1:nimages
        R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
        theta_opt = atan2d(R_tmp(2,1), R_tmp(1,1));
        image_set_aligned(:,:,i) = imrotate(image_set_subsampled(:,:,i), -theta_opt, 'crop');
    end
    
    %
    if corr(time_subsampled, embed_coord(:,1), 'type','spearman') < 0
        embed_coord(:,1) = -embed_coord(:,1);
    end
    
    [~, I] = sort(embed_coord(:,1));
    % figure;
    % for i=1:nimages
    %     subplot(6, 7,i)
    %     imshow(image_set_aligned(:,:,I(i)))
    % end
    
    % make_fig(4,4);
    % plot(compute_ranks(time_subsampled), compute_ranks(embed_coord(:,1)),'.')
    % title(sprintf('corr = %2.2f', corr(compute_ranks(time_subsampled), compute_ranks(embed_coord(:,1)))))
    % xlabel('true rank')
    % ylabel('recovered rank')
    % set(gca, 'xtick', [0 20 40])
    % set(gca, 'ytick', [0 20 40])
    % axis([0 50 0 50])
    % axis square
    
    rank_corr(j) = corr(compute_ranks(time_subsampled), compute_ranks(embed_coord(:,1)));
    
%     if rank_corr(j) > 0.80 && rank_corr(j) < 0.81
%         make_fig(4,4);
%         plot(compute_ranks(time_subsampled), compute_ranks(embed_coord(:,1)),'.')
%         xlabel('true rank')
%         ylabel('recovered rank')
%         set(gca, 'xtick', [0 20 40])
%         set(gca, 'ytick', [0 20 40])
%         axis([0 50 0 50])
%         axis square
%         saveas(gcf, 'shuffled_movie_rank_corr.pdf');
%         
%         nprint_images = 10;
%         
%         make_fig(17, 17/nprint_images);
%         [~, I] = sort(embed_coord(:, 1));
%         for i1=1:nprint_images
%             i = I(round(nimages*i1/nprint_images));
%             im_tmp = image_set_aligned(:,:,i);
%             make_subplot(nprint_images, 1, 0.01, i1);
%             imshow(imrotate(im_tmp, 140, 'crop'));
%             text(npixels/2, npixels/2, sprintf('%2.1f min \n movie=%d', time_subsampled(i), idx(i)),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
%             
%         end
%         saveas(gcf, 'shuffled_movie_registered_ordered.pdf');
%         
%         return
%     end
end

figure;
hist(rank_corr, 20)
title(sprintf('median rank corr = %2.2f', median(rank_corr)))



