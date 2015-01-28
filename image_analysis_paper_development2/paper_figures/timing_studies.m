clear all
close all

%%

% addpath ../code_for_distribution

%%
image_name = 'image';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
% npixels = 100;
dim = 2;

channel_weight = 1;
channel_blur = 0;
channel_normalize = 0;
channel_mean_center = 1;
resize_image = false;

ang_dis = 10;

eps_scale = 0.25;

image_dir = '../code_for_distribution/zebrafish';

npixels_zebrafish = [50 100 150 200 300 400];
cpu_times_npixels_zebrafish = zeros(size(npixels_zebrafish));

for i=1:length(npixels_zebrafish)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels_zebrafish(i), dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, ang_dis);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_npixels_zebrafish(i) = toc;
end


%%

npixels = 100;
nimages_zebrafish = [20 40 60 80 100 120];
cpu_times_nimages_zebrafish = zeros(size(nimages_zebrafish));
for i=1:length(nimages_zebrafish)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages_zebrafish(i), nstack, npixels, dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, ang_dis);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_nimages_zebrafish(i) = toc;
end

%%
npixels = 100;
nimages = 120;
nrot_zebrafish = [12 24 36 72];
cpu_times_nrot_zebrafish = zeros(size(nrot_zebrafish));
for i=1:length(nrot_zebrafish)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, 360/nrot_zebrafish(i));
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_nrot_zebrafish(i) = toc;
end

%%
image_name = 'emb';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
dim = 2;

channel_weight = [0.5 1 1];
channel_blur = [0.05 0.05 0.05];
channel_normalize = [1 0 0];
channel_mean_center = [1 0 0];
resize_image = true;

ang_dis = 10;

eps_scale = 0.25;

image_dir = '../code_for_distribution/drosophila_fixed_images';

npixels_drosophila = [50 100 150 200 300 400];
cpu_times_npixels_drosophila = zeros(size(npixels_drosophila));

for i=1:length(npixels_drosophila)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels_drosophila(i), dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, ang_dis);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_npixels_drosophila(i) = toc;
end


%%

npixels = 100;
nimages_drosophila = [20 40 60 80 100 120];
cpu_times_nimages_drosophila = zeros(size(nimages_drosophila));

for i=1:length(nimages_drosophila)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages_drosophila(i), nstack, npixels, dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, ang_dis);
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_nimages_drosophila(i) = toc;
end

%%

npixels = 100;
nimages = 120;
nrot_drosophila = [12 24 36 72];
cpu_times_nrot_drosophila = zeros(size(nrot_drosophila));

for i=1:length(nrot_drosophila)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, 360/nrot_drosophila(i));
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_nrot_drosophila(i) = toc;
end

%%

fontsize = 8;
figsize = 5;

make_fig(figsize, figsize);
loglog(npixels_zebrafish, cpu_times_npixels_zebrafish, '.b', 'markersize', 10);
hold on
loglog(npixels_drosophila, cpu_times_npixels_drosophila, '.r', 'markersize', 10);
xlabel('number of pixels', 'fontsize', fontsize)
ylabel('CPU time (seconds)', 'fontsize', fontsize)
h = legend('grayscale','color','location','southeast');
set(h, 'fontsize', 4);
set(gca, 'fontsize', fontsize)
set(gca, 'xtick', [100 1000])
set(gca, 'xticklabel', {'100x100'; '1000x1000'});
x_data = linspace(50, 1500, 10);
P = polyfit(log(npixels_zebrafish(2:end)), log(cpu_times_npixels_zebrafish(2:end)),1);
P(1)
loglog(x_data, exp(polyval(P, log(x_data))), '--b');

P = polyfit(log(npixels_drosophila(2:end)), log(cpu_times_npixels_drosophila(2:end)),1);
P(1)
loglog(x_data, exp(polyval(P, log(x_data))), '--r');
axis tight
axis square
saveas(gcf, 'timing_studies_npixels.pdf');

make_fig(figsize, figsize);
loglog(nimages_zebrafish, cpu_times_nimages_zebrafish, '.b', 'markersize', 10);
hold on
loglog(nimages_drosophila, cpu_times_nimages_drosophila, '.r', 'markersize', 10);
xlabel('number of images', 'fontsize', fontsize)
ylabel('CPU time (seconds)', 'fontsize', fontsize)
h = legend('grayscale','color','location','southeast');
set(h, 'fontsize', 4);
set(gca, 'fontsize', fontsize)

x_data = linspace(20, 1000, 10);
P = polyfit(log(nimages_zebrafish(2:end)), log(cpu_times_nimages_zebrafish(2:end)),1);
P(1)
loglog(x_data, exp(polyval(P, log(x_data))), '--b');

P = polyfit(log(nimages_drosophila(2:end)), log(cpu_times_nimages_drosophila(2:end)),1);
P(1)
loglog(x_data, exp(polyval(P, log(x_data))), '--r');

axis tight
axis square
saveas(gcf, 'timing_studies_nimages.pdf');


make_fig(figsize, figsize);
loglog(360./nrot_zebrafish, cpu_times_nrot_zebrafish, '.b', 'markersize', 10);
hold on
loglog(360./nrot_drosophila, cpu_times_nrot_drosophila, '.r', 'markersize', 10);
xlabel('angular discretization', 'fontsize', fontsize)
ylabel('CPU time (seconds)', 'fontsize', fontsize)
h = legend('grayscale','color','location','southwest');
set(h, 'fontsize', 4);
set(gca, 'fontsize', fontsize)
set(gca, 'xtick', [1 10])
set(gca, 'xticklabel', {'1°'; '10°'})

x_data = linspace(12, 360, 10);
P = polyfit(log(nrot_zebrafish(2:end)), log(cpu_times_nrot_zebrafish(2:end)),1);
P(1)
loglog(360./x_data, exp(polyval(P, log(x_data))), '--b');

P = polyfit(log(nrot_drosophila(2:end)), log(cpu_times_nrot_drosophila(2:end)),1);
P(1)
loglog(360./x_data, exp(polyval(P, log(x_data))), '--r');
axis tight
axis square
saveas(gcf, 'timing_studies_nrot.pdf');

