clear all
close all

%%

% addpath ../code_for_distribution

%%
image_name = 'image';
image_ext = 'tif';
stack_name = '';
nimages = 40;
nstack = 0;
npixels = 100;
dim = 2;

channel_weight = 1;
channel_blur = 0.05;
channel_normalize = 1;
channel_mean_center = 0;
resize_image = false;

ang_dis = 10;

eps_scale = 0.25;

nexpt = 5;
rank_corr = zeros(nexpt, 1);
theta_err = zeros(nexpt, 1);

for i=1:nexpt
    image_dir = sprintf('../code_for_distribution/drosophila_live_imaging/expt%02d', i);
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
    %         plot_images(images_raw, dim)
    
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    %         plot_images(images, dim)
    %     figure;
    %     for j=1:nimages
    %         subplot(4, 10, j)
    %         imshow(edge(images(:,:,j), 'canny'));
    %     end
    %     continue
    
    [R, W] = compute_pairwise_alignments(images, ang_dis);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    %         plot_images(images_analyzed, dim)
    
    fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
    time = fscanf(fid, '%f');
    fclose(fid);
    
    %     figure;
    %     plot(tiedrank(time), tiedrank(embed_coord),'.')
    
    rank_corr(i) = corr(time, embed_coord,'type','spearman');
    
    fid = fopen(sprintf('%s/angles.txt', image_dir), 'r');
    angles = fscanf(fid, '%f');
    fclose(fid);
    
    angles_opt = R_to_theta(R_opt);
    [~, angles_opt] = angle_corr(angles, angles_opt);
    
    %     figure;
    %     plot(angles, angles_opt, '.')
    
    %     corr(angles, angles_opt)
    
    theta_err(i) = mean(abs(angles - angles_opt));
    
    if i == 1
        nsubimages = 12;
        fontsize = 8;
        if corr(time, embed_coord) < 0
            embed_coord = -embed_coord;
        end
        [~, I] = sort(embed_coord(1:nsubimages));
        make_fig(17, 17/nsubimages);
        for j=1:nsubimages
            make_subplot(nsubimages, 1, 0.01, j);
            imshow(images(:,:,j));
        end
        saveas(gcf, 'drosophila_live_imaging_scrambled.pdf');
        
        make_fig(17, 17/nsubimages);
        for j=1:nsubimages
            make_subplot(nsubimages, 1, 0.01, j);
            imshow(imrotate(images_registered(:,:,I(j)), 50, 'crop'));
            %             text(npixels/2, npixels/2, sprintf('%2.1f min', time(I(j))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
            text(npixels, npixels, sprintf('%2.1f min', time(I(j))),'color',0.95*ones(1,3),'HorizontalAlignment','right','VerticalAlignment','bottom', 'fontsize', 4)
        end
        saveas(gcf, 'drosophila_live_imaging_ordered.pdf');
        
        make_fig(5, 5);
        plot(angles, angles_opt,'.k')
        xlabel('true angle', 'fontsize', fontsize)
        ylabel('recovered angle', 'fontsize', fontsize)
        set(gca, 'fontsize', fontsize)
        set(gca, 'xtick', [0 180 360]);
        set(gca, 'ytick', [0 180 360]);
        set(gca, 'xticklabel', {'0°'; '180°'; '360°'});
        set(gca, 'yticklabel', {'0°'; '180°'; '360°'});
        axis([0 400 0 400])
        axis square
        saveas(gcf, 'drosophila_live_imaging_angle_corr.pdf');
        
        make_fig(5, 5);
        plot(tiedrank(time), tiedrank(embed_coord),'.k')
        xlabel('true rank', 'fontsize', fontsize)
        ylabel('recovered rank', 'fontsize', fontsize)
        set(gca, 'fontsize', fontsize)
        set(gca, 'xtick', [0 20 40])
        set(gca, 'ytick', [0 20 40])
        set(gca, 'xticklabel', {'0 '; '  20 '; '  40 '});
        set(gca, 'yticklabel', {'0 '; '  20 '; '  40 '});
        axis([0 50 0 50])
        axis square
        saveas(gcf, 'drosophila_live_imaging_rank_corr.pdf');
        
        ncomps = 20;
        [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
        
        make_fig(8,8)
        plot(abs(D2), '.')
        xlabel('embedding coordinate')
        ylabel('product of eigenvalues')
        axis square
        saveas(gcf, 'drosophila_live_eval_spectrum.pdf');
        
        for k = 2:4
            make_fig(3,3)
            plot(embed_coord(:,1), embed_coord(:,k), '.')
            xlabel('embed coord 1')
            ylabel(sprintf('embed coord %d', k))
            axis tight
            axis square
            set(gca, 'xtick', [])
            set(gca, 'ytick', [])
            saveas(gcf, sprintf('drosophila_live_evec_corr%d.pdf', k-1));
        end
        return
        
        
        
        
    end
end

rank_corr
theta_err

%%

channel_weight = 1;
channel_blur = 0.05;
channel_normalize = 1;
channel_mean_center = 1;
resize_image = true;

eps_scale = 0.5;

% images_raw_all = zeros(npixels, npixels, nimages, nexpt);
images_all = zeros(npixels, npixels, nimages, nexpt);

for i=1:nexpt
    image_dir = sprintf('../code_for_distribution/drosophila_live_imaging/expt%02d', i);
    
    [images_raw_tmp, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
    
    fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
    time = fscanf(fid, '%f');
    fclose(fid);
    
    [~, idx] = sort(time);
    images_raw_tmp = images_raw_tmp(:,:,idx);
    
    images_tmp = apply_image_functions(images_raw_tmp, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    images_all(:,:,:,i) = images_tmp;
    
end

time = (1:nimages)';

nreps = 50;
rank_corr = zeros(nreps, 1);
for k=1:nreps
    rng(21*k)
    movie_idx = randi(nexpt, [nimages 1]);
    
    for i=1:nexpt
        I = find(movie_idx == i);
        images(:,:,I) = images_all(:,:,I,i);
    end
    
    [R, W] = compute_pairwise_alignments(images, nrot);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    %     images_registered = register_all_images(images, R_opt);
    %
    %     images_analyzed = order_all_images(images_registered, embed_coord);
    %
    %     plot_images(images_analyzed, dim)
    %
    %     figure;
    %     plot(tiedrank(time), tiedrank(embed_coord),'.')
    
    rank_corr(k) = abs(corr(time, embed_coord,'type','spearman'));
end

median(rank_corr)


%% bootstrap
% rng(321);
% nsamples = 10;
% nreps = 50;
% sample_vec = round(linspace(10, nimages, nsamples));
% theta_err = zeros(nsamples, nreps);
% rank_corr = zeros(nsamples, nreps);
%
% i = 1;
% image_dir = sprintf('../code_for_distribution/drosophila_live_imaging/expt%02d', i);
% [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
%
% images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
%
% [R, W] = compute_pairwise_alignments(images, nrot);
%
% fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
% time = fscanf(fid, '%f');
% fclose(fid);
%
% fid = fopen(sprintf('%s/angles.txt', image_dir), 'r');
% angles = fscanf(fid, '%f');
% fclose(fid);
%
% for j = 1:nsamples
%
%     j
%
%     for j2 = 1:nreps
%
%         sample_idx = randi(nimages, sample_vec(j), 1);
%
%         R_idx = reshape([2*sample_idx'-1; 2*sample_idx'], [], 1);
%
%         %         [R, W] = compute_pairwise_alignments(images(:,:,sample_idx), nrot);
%
%         % compute optimal rotations + embedding coordinates using vector diffusion maps
%         [R_opt, embed_coord, D2] = vdm(R(R_idx, R_idx), W(sample_idx, sample_idx), eps_scale, ncomps);
%
%         % register images using optimal rotations
%         %         images_registered = register_all_images(images(:,:,sample_idx), R_opt);
%
%         %         images_analyzed = order_all_images(images_registered, embed_coord);
%
%
%         rank_corr(j, j2) = abs(corr(time(sample_idx), embed_coord,'type','spearman'));
%
%         angles_opt = R_to_theta(R_opt);
%         [theta_err(j, j2), ~] = angle_corr(angles(sample_idx), angles_opt);
%
%     end
% end
%
% figure;
% plot(sample_vec, median(theta_err, 2), '.')
% set(gca, 'ylim', [0 1])
%
% figure;
% plot(sample_vec, median(rank_corr, 2), '.')
% set(gca, 'ylim', [0 1])