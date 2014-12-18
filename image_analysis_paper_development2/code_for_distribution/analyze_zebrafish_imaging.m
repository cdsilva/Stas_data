clear all
close all

%% read in images

% number of pixels
npixels = 100;
% time step for live imaging
dt = 1;

% open live imaging movie file
obj = VideoReader('zebrafish/Zebrafish.mp4');
% read in images
images_raw = obj.read;

% get number of images
nimages = get(obj, 'NumberOfFrames');

% store times of images
time = (1:nimages)' * dt;

% process images
images = zeros(npixels, npixels, nimages, 'uint8');
for i=1:nimages
    
    % extract relevant channel; resize to number of pixels
    im_tmp = imresize(images_raw(:,:, 1, i), [npixels npixels]);
    
    % crop image
%     im_tmp = crop_image(im_tmp, nuclear_channel);
    
    % normalize signal
%     im_tmp = normalize_nuclear_signal(im_tmp);

    images(:, :, i) = im_tmp;
end

%%

idx = 53:172;
images = images(:,:,idx);
time = time(idx);
nimages = length(idx);

make_subplot = @(j) subplot(8, 15, j);

%% show raw images

figure; 
for j=1:nimages
    make_subplot(j)
    imshow(images(:, :, j));
end

%% rotate and shuffle images

% initialize random number generator
rng(213);

% random angles to rotate images
theta = 360*rand(nimages, 1);
% indices to shuffle images
shuffle_idx = randperm(nimages);

for i=1:nimages
    images(:, :, i) = imrotate(images(:, :, i), theta(i), 'crop');
end

images = images(:, :, shuffle_idx);
time = time(shuffle_idx);
theta = theta(shuffle_idx);

%% show rotated and shuffled images

figure;
for j=1:nimages
    make_subplot(j)
    imshow(images(:, :, j));
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

if corr(embed_coord(:, 1), (1:nimages)') < 0
    embed_coord(:, 1) = -embed_coord(:, 1);
end

[~, I] = sort(embed_coord(:, 1));
images_registered_ordered = images_registered(:,:,I);

%% compare angles

figure;
plot(mod(theta, 360), mod(R_to_theta(R_opt), 360), '.k')
xlabel('true angle')
ylabel('recovered angle')

%% compare ranks

figure;
plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.k')
xlabel('true rank')
ylabel('recovered rank')

%% show registered and ordered images

figure;
for j=1:nimages
    make_subplot(j)
    imshow(images_registered_ordered(:, :, j));
end
