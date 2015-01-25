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
channel_normalize = 0;
channel_mean_center = 1;
resize_image = false;

nrot = 36;

eps_scale = 0.25;

for i=1:5
    image_dir = sprintf('../code_for_distribution/drosophila_live_imaging/expt%02d', i);
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
%     plot_images(images_raw, dim)
    
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
%     plot_images(images, dim)
    
    [R, W] = compute_pairwise_alignments(images, nrot);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
%     plot_images(images_analyzed, dim)
    
    fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
    time = fscanf(fid, '%f');
    fclose(fid);
    
    figure;
    plot(tiedrank(time), tiedrank(embed_coord),'.')
    
    corr(time, embed_coord,'type','spearman')
    
    fid = fopen(sprintf('%s/angles.txt', image_dir), 'r');
    angles = fscanf(fid, '%f');
    fclose(fid);
    
    angles_opt = R_to_theta(R_opt);
    angles_opt = angles_opt + mean(atan2d(sind(angles-angles_opt), cosd(angles-angles_opt)));
    angles_opt = mod(angles_opt, 360);
    idx = find(angles - angles_opt > 180);
    angles_opt(idx) = angles_opt(idx) + 360;
    idx = find(angles - angles_opt < -180);
    angles_opt(idx) = angles_opt(idx) - 360;
    
    figure;
    plot(angles, angles_opt, '.')
    
    corr(angles, angles_opt)
    
end