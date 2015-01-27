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

nrot = 36;

eps_scale = 0.25;

image_dir = '../code_for_distribution/zebrafish';

npixels_zebrafish = [50 100 150 200 300 400];
cpu_times_npixels_zebrafish = zeros(size(npixels_zebrafish));

for i=1:length(npixels_zebrafish)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels_zebrafish(i), dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, nrot);
    
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
    
    [R, W] = compute_pairwise_alignments(images, nrot);
    
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
    
    [R, W] = compute_pairwise_alignments(images, nrot_zebrafish(i));
    
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

nrot = 36;

eps_scale = 0.25;

image_dir = '../code_for_distribution/drosophila_fixed_images';

npixels_drosophila = [50 100 150 200 300 400];
cpu_times_npixels_drosophila = zeros(size(npixels_drosophila));

for i=1:length(npixels_drosophila)
    
    [images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels_drosophila(i), dim);
    
    tic
    images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
    
    [R, W] = compute_pairwise_alignments(images, nrot);
    
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
    
    [R, W] = compute_pairwise_alignments(images, nrot);
    
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
    
    [R, W] = compute_pairwise_alignments(images, nrot_drosophila(i));
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
    
    % register images using optimal rotations
    images_registered = register_all_images(images, R_opt);
    
    images_analyzed = order_all_images(images_registered, embed_coord);
    
    cpu_times_nrot_drosophila(i) = toc;
end

%%

make_fig(8,8);
loglog(npixels_zebrafish, cpu_times_npixels_zebrafish, '.-b');
hold on
loglog(npixels_drosophila, cpu_times_npixels_drosophila, '.-r');
xlabel('number of pixels')
ylabel('CPU time')
legend('grayscale','color','location','northwest')
set(gca, 'xtick', [50 100 200 400])
set(gca, 'ytick', [10 60 360])
set(gca, 'YTickLabel', {'10 sec'; '1 min'; '6 min'})
axis tight
axis square
saveas(gcf, 'timing_studies_npixels.pdf');

make_fig(8,8);
loglog(nimages_zebrafish, cpu_times_nimages_zebrafish, '.-b');
hold on
loglog(nimages_drosophila, cpu_times_nimages_drosophila, '.-r');
xlabel('number of images')
ylabel('CPU time')
legend('grayscale','color','location','northwest')
set(gca, 'xtick', [20 40 80 160])
set(gca, 'ytick', [2 5 10 20])
set(gca, 'YTickLabel', {'2 sec'; '5 sec'; '10 sec'; '20 sec'})
axis tight
axis square
saveas(gcf, 'timing_studies_nimages.pdf');

make_fig(8,8);
loglog(nrot_zebrafish, cpu_times_nrot_zebrafish, '.-b');
hold on
loglog(nrot_drosophila, cpu_times_nrot_drosophila, '.-r');
xlabel('number of rotations')
ylabel('CPU time')
legend('grayscale','color','location','northwest')
set(gca, 'xtick', [10 20 40 80])
set(gca, 'ytick', [10 20 40])
set(gca, 'YTickLabel', {'10 sec'; '20 sec'; '40 sec'})
axis tight
set(gca, 'xlim', [10 80])
axis square
saveas(gcf, 'timing_studies_nrot.pdf');

