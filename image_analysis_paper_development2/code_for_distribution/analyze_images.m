clear all
close all

%% read in images

% number of pixels
npixels = 100;
% number of images
nimages = 120;
% image directory
image_dir = 'drosophila_fixed_images';
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
    
    % normalize+scale nuclear signal
    im_tmp(:,:,nuclear_channel) = normalize_nuclear_signal(im_tmp(:,:,nuclear_channel));
    
    % crop image
    im_tmp = crop_image(im_tmp, nuclear_channel);
    
    im_tmp(:, :, nuclear_channel) = immultiply(im_tmp(:, :, nuclear_channel), 0.5);
    
    images(:, :, :, i) = im_tmp;
end

% show raw images
figure;
for i=1:nimages
    subplot(8, 15, i)
    imshow(images(:,:,:,i))
end


% compute pairwise alignments + distances
nrot = 36;
show_waitbar = true;
[R, W] = compute_pairwise_alignments(images, nrot, show_waitbar);

% compute optimal rotations + embedding coordinates using vector diffusion maps
eps_scale = sqrt(0.1);
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);
[~, I] = sort(embed_coord(:, 1));
images_registered_ordered = images_registered(:,:,:,I);


% show registered and ordered images
figure;
for i=1:nimages
    subplot(8, 15, i)
    imshow(images_registered_ordered(:,:,:,i));
end

figure;
plot(compute_ranks(embed_coord(:,1)),'.')
xlabel('expert rank')
ylabel('algorithm rank')

%%
% number of pixels
npixels = 100;
% number of images
nimages = 40;
% image directory
image_dir = 'drosophila_live_imaging/expt05';
% image prefix
image_name = 'image';
% number of channels in images
nchannels = 1;

% process images
images = zeros(npixels, npixels, nimages, 'uint8');
for i=1:nimages
    
    % read image
    im_tmp = imread(sprintf('%s/%s%02d.tif', image_dir, image_name, i));
    
    % extract relevant channel; resize to number of pixels
    im_tmp = imresize(im_tmp, [npixels npixels]);
    
    
    % normalize+scale nuclear signal
    im_tmp = normalize_nuclear_signal(im_tmp);
    
    % crop image
    im_tmp = crop_image(im_tmp);
    
    images(:, :, i) = im_tmp;
end

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

% show raw images
figure;
for i=1:nimages
    subplot(4, 10, i)
    imshow(images(:,:,i))
end

% compute pairwise alignments + distances
nrot = 36;
show_waitbar = true;
[R, W] = compute_pairwise_alignments(images, nrot, show_waitbar);

% compute optimal rotations + embedding coordinates using vector diffusion maps
eps_scale = 0.25;
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);
[~, I] = sort(embed_coord(:, 1));
images_registered_ordered = images_registered(:,:,I);


% show registered and ordered images
figure;
for i=1:nimages
    subplot(4, 10, i)
    imshow(images_registered_ordered(:,:,i));
end

figure;
plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.')
xlabel('expert rank')
ylabel('algorithm rank')

%%

% number of pixels
npixels = 100;
% number of images
nimages = 120;
% image directory
image_dir = 'zebrafish';
% image prefix
image_name = 'image';
% number of channels in images
nchannels = 1;

% process images
images = zeros(npixels, npixels, nimages, 'uint8');
for i=1:nimages
    
    % read image
    im_tmp = imread(sprintf('%s/%s%02d.tif', image_dir, image_name, i));
    
    % extract relevant channel; resize to number of pixels
    im_tmp = imresize(im_tmp, [npixels npixels]);
    
    % crop image
    im_tmp = crop_image(im_tmp);
    
    % normalize+scale nuclear signal
    im_tmp = normalize_nuclear_signal(im_tmp);
    
    images(:, :, i) = im_tmp;
end

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

% show raw images
figure;
for i=1:nimages
    subplot(12, 10, i)
    imshow(images(:,:,i))
end

% compute pairwise alignments + distances
nrot = 36;
show_waitbar = true;
[R, W] = compute_pairwise_alignments(images, nrot, show_waitbar);

% compute optimal rotations + embedding coordinates using vector diffusion maps
eps_scale = 0.25;
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);
[~, I] = sort(embed_coord(:, 1));
images_registered_ordered = images_registered(:,:,I);


% show registered and ordered images
figure;
for i=1:nimages
    subplot(12, 10, i)
    imshow(images_registered_ordered(:,:,i));
end

figure;
plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.')
xlabel('expert rank')
ylabel('algorithm rank')

%%

% number of pixels
npixels = 100;

nchannels = 3;

nstacks = 46;
nimages_stack = 21;

% image directory
image_dir = 'wing_disc_stacks';
stack_dir_prefix = 'disc';
image_prefix = 'z';

image_stacks = zeros(npixels, npixels, nchannels, nimages_stack, nstacks, 'uint8');


for i=1:nstacks
    for j=1:nimages_stack
        
        im_tmp = imread(sprintf('%s/%s%02d/%s%02d.tif', image_dir, stack_dir_prefix, i, image_prefix, j));
        im_tmp = imresize(im_tmp, [npixels npixels]);
        
        image_stacks(:,:,:,j,i) = im_tmp;
        
    end
    image_stacks(:,:,:,:,i) = crop_zstack(image_stacks(:,:,:,:,i));
    
end

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

figure;
for i=1:nstacks
    subplot(6, 8, i)
    imshow(max(image_stacks(:,:,:,:,i), [], 4));
end

% compute pairwise alignments + distances
nrot = 36;
show_waitbar = true;
[R, W] = compute_pairwise_alignments(squeeze(max(image_stacks, [], 4)), nrot, show_waitbar);

% compute optimal rotations + embedding coordinates using vector diffusion maps
eps_scale = 0.25;
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);
images_registered = register_all_images(image_stacks, R_opt);


figure;
for i=1:nstacks
    subplot(6, 8, i)
    imshow(immultiply(max(images_registered(:,:,:,:,i), [], 4), 1.5));
end

W = squareform(pdist(reshape(double(images_registered), [], nstacks)')).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

[~, I] = sort(V(:,2));
figure;
for i=1:nstacks
    subplot(6, 8, i)
    imshow(immultiply(max(images_registered(:,:,:,:,I(i)), [], 4), 1.5));
end

figure; 
plot(time(I),'.')