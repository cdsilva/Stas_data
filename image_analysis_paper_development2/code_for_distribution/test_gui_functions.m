clear all
close all

%%
image_dir = 'drosophila_fixed_images';
image_name = 'emb';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
npixels = 100;
dim = 2;

channel_weight = [0.5 1 1];
channel_blur = [0.05 0 0];
channel_normalize = [1 0 0];
channel_mean_center = [1 0 0];
resize_image = true;

nrot = 36;

eps_scale = 0.25;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(images, nrot);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

images_analyzed = order_all_images(images_registered, embed_coord);

plot_images(images_analyzed, dim)

time = (1:nimages)';
figure;
plot(tiedrank(time), tiedrank(embed_coord),'.')

corr(time, embed_coord,'type','spearman')

%%
image_dir = 'drosophila_live_imaging/expt01';
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

nrot = 36;

eps_scale = 0.2;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(images, nrot);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

images_analyzed = order_all_images(images_registered, embed_coord);

plot_images(images_analyzed, dim)

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

figure;
plot(tiedrank(time), tiedrank(embed_coord),'.')

corr(time, embed_coord,'type','spearman')

%%

image_dir = 'wing_disc_stacks';
image_name = 'disc';
image_ext = 'tif';
stack_name = 'z';
nimages = 46;
nstack = 21;
npixels = 100;
dim = 3;

channel_weight = [1 1 1];
channel_blur = [0 0 0];
channel_normalize = [0 0 0];
channel_mean_center = [0 0 1];
resize_image = false;

nrot = 36;

eps_scale = 0.5;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(squeeze(max(images, [], ndims(images)-1)), nrot);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

plot_images(images_registered, dim)

W = squareform(pdist(reshape(double(images_registered), [], nimages)')).^2;
[embed_coord, D2] = dm(W, eps_scale, ncomps);

images_analyzed = order_all_images(images_registered, embed_coord);

plot_images(immultiply(images_analyzed, 1.5), dim)

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

figure;
plot(tiedrank(time), tiedrank(embed_coord),'.')

corr(time, embed_coord,'type','spearman')

% image_dir2 = 'wing_disc_stacks2';
% save_images(images_analyzed, dim, image_dir2, image_name, image_ext, stack_name)

%%

image_dir = 'zebrafish';
image_name = 'image';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
npixels = 100;
dim = 2;

channel_weight = 1;
channel_blur = 0.05;
channel_normalize = 1;
channel_mean_center = 1;
resize_image = false;

nrot = 36;

eps_scale = 0.2;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(images, nrot);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

images_analyzed = order_all_images(images_registered, embed_coord);

plot_images(images_analyzed, dim)

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

figure;
plot(tiedrank(time), tiedrank(embed_coord),'.')

corr(time, embed_coord,'type','spearman')
