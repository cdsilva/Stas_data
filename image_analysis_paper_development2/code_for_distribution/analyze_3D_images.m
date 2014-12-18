clear all
close all

%% read in images

% number of pixels
npixels = 100;
npixels2 = 300;
% number of images
nimages = 45;
% image directory
image_dir = '../../data_3d/3D_data';
% image prefix
image_name = 'emb';
% channel where nuclear signal is stored (1=red, 2=green, 3=blue)
nuclear_channel = 1;
% number of channels in images
nchannels = 3;

% process images
images = zeros(npixels, npixels, npixels2, nchannels, nimages, 'uint8');
for i=1:nimages
    
    i 
    j = 1;
    images_tmp = [];
    
    while exist(sprintf('%s/%s%02d/%s%02d_z%03d.tif', image_dir, image_name, i, image_name, i, j), 'file')
        
        % read image
        im_tmp = imread(sprintf('%s/%s%02d/%s%02d_z%03d.tif', image_dir, image_name, i, image_name, i, j));
        
        % extract relevant channel; resize to number of pixels
        im_tmp = imresize(im_tmp, [npixels2 npixels]);
        
        images_tmp = cat(4, images_tmp, im_tmp);
        
        j = j + 1;
    end
    
    images_tmp = permute(images_tmp, [2 4 1 3]);
    images_tmp = imresize(images_tmp, [npixels npixels]);
    
    images(:, :, :, :, i) = images_tmp;
end


%%

% [X, Y, Z] = meshgrid((1:npixels)-(mean(1:npixels)), (1:npixels)-(mean(1:npixels)), (1:npixels2)-(mean(1:npixels2)));

% I = find((X/(npixels/2)).^2+(Y/(npixels/2)).^2+(Z/(npixels2/2)).^2 < 0.7);

% images_tmp = images;
% for i=1:nimages
%     for j=1:3
%         im_tmp = images_tmp(:,:,:,j,i);
%         im_tmp(I) = 0;
%         images_tmp(:,:,:,j,i) = im_tmp;
%     end
% end

        
% images_tmp(25:75, 25:75, 75:225, :, :) = 0;

% figure;
% for i=1:nimages
%     subplot(9, 5, i)
%     imshow(squeeze(max(images_tmp(:,:,:,1,i), [], 2)));
% end

figure;
for i=1:nimages
    subplot(9, 5, i)
    imshow(squeeze(images(round(npixels/2),:,:,1,i)));
end

figure;
for i=1:nimages
    subplot(9, 5, i)
    imshow(immultiply(squeeze(max(images(:,:,:,1,i), [], 1)),2));
end

return


%% compute pairwise alignments + distances

% number of angles used to compute the pairwise alignments
nrot = 12;
% do not consider shifts
nshifts = 0;
shift_max = 0;
% show progress bar
show_waitbar = true;

images_tmp = images;
% images_tmp(:,:,:,2,:) = 0;
% images_tmp(:,:,:,3,:) = 0;
% images_tmp(:,:,:,1,:) = immultiply(images_tmp(:,:,:,1,:), 2);

[R, W] = compute_pairwise_alignments(images_tmp, nrot, shift_max, nshifts, show_waitbar);

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
images_registered_ordered = images_registered(:,:,:,:,I);

%% show registered and ordered images

figure;
for i=1:nimages
    subplot(9, 5, i)
    imshow(immultiply(squeeze(max(images(:,:,:,1,i), [], 1)),2));
end

return
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
