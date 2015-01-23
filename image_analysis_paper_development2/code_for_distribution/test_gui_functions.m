clear all
close all

image_dir = 'drosophila_fixed_images';
image_name = 'emb';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
npixels = 100;
dim = 2;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);

subplot_dim1 = 12;
subplot_dim2 = 10;

plot_images(images_raw, dim, subplot_dim1, subplot_dim2)

channel_weight = [0.5 1 1];
channel_blur = [0.05 0 0];
channel_normalize = [1 0 0];
channel_mean_center = [1 0 0];
resize_image = true;

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);

plot_images(images, dim, subplot_dim1, subplot_dim2)

nrot = 36;
show_waitbar = true;
[R, W] = compute_pairwise_alignments(images, nrot, show_waitbar);


% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
eps_scale = 0.3;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

[~, I] = sort(embed_coord);

if ndims(images_registered) == 3
    images_analyzed = images_registered(:,:,I);
elseif ndims(images_registered) == 4
    images_analyzed = images_registered(:,:,:,I);
elseif ndims(images_registered) == 5
    images_analyzed = images_registered(:,:,:,I);
end

plot_images(images_analyzed, dim, subplot_dim1, subplot_dim2)

corr((1:nimages)', embed_coord,'type','spearman')

