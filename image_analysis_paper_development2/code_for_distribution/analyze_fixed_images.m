clear all
close all

%% read in images

% number of pixels
npixels = 100;
% number of images
nimages = 120;
% image directory
image_dir = 'fixed_images';
% image prefix 
image_name = 'emb';
% channel where nuclear signal is stored (1=red, 2=green, 3=blue)
nuclear_channel = 1;
% number of channels in images
nchannels = 3;

% process images
images = zeros(npixels, npixels, nchannels, nimages, 'uint8');
for i=1:nimages
    
    % read image
    im_tmp = imread(sprintf('%s/%s%02d.tif', image_dir, image_name, i));
    
    % extract relevant channel; resize to number of pixels
    im_tmp = imresize(im_tmp, [npixels npixels]);
    
    % crop image
    im_tmp = crop_image(im_tmp, nuclear_channel);
    
    % normalize+scale nuclear signal
    im_tmp(:,:,nuclear_channel) = normalize_nuclear_signal(im_tmp(:,:,nuclear_channel));
    im_tmp(:, :, nuclear_channel) = immultiply(im_tmp(:, :, nuclear_channel), 0.5);

    images(:, :, :, i) = im_tmp;
end

%% show raw images

figure;
for i=1:nimages
    subplot(8, 15, i)
    imshow(images(:,:,:,i))
end

%% compute pairwise alignments + distances

% number of angles used to compute the pairwise alignments
nrot = 36;
% do not consider shifts
nshifts = 0;
shift_max = 0;
% show progress bar
show_waitbar = true;

[R, W] = compute_pairwise_alignments(images, nrot, shift_max, nshifts, show_waitbar);

%% compute optimal rotations + embedding coordinates using vector diffusion maps

% use default kernel scale
eps = 0;
% only compute first embedding coordinate
ncomps = 1;

[R_opt, embed_coord, D2] = vdm(R, W, eps, ncomps);

%% register images using optimal rotations

images_registered = register_all_images(images, R_opt);

%% order images using first embedding coordinate

[~, I] = sort(embed_coord(:, 1));
images_registered_ordered = images_registered(:,:,:,I);


%% show registered and ordered images

figure;
for i=1:nimages
    subplot(8, 15, i)
    imshow(images_registered_ordered(:,:,:,i));
end

%% compare ranks
figure;
plot(compute_ranks(embed_coord(:,1)),'.')
xlabel('expert rank')
ylabel('algorithm rank')

%% calculate average trajectory

% number of images in average trajectory
nstages = 10;
% kernel scale used for computing averages
window_scale = 3;

avg_images = compute_average_trajectory(images_registered_ordered, nstages, window_scale);

%% show average trajectory

figure;
for i=1:nstages
    subplot(1, nstages, i);
    imshow(avg_images(:,:,:,i))
end
