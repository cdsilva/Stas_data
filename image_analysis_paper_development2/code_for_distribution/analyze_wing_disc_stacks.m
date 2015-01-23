clear all
close all

%% read in images

% number of pixels
npixels = 100;

nchannels = 3;

nimages_stack_UD = 10;

% image directory
image_dir = 'wing_disc_stacks';

folders = {'72-73H'; '76,5-77,5H'; '79-80H'; '89-90H'; '100-101H'; '110,5-111,5H'};
time_for_folders = [72 76.5 79 89 100 110.5];

image_stacks = [];
time = [];

for i=1:length(folders)
    nstacks = length(dir(sprintf('%s/%s', image_dir, folders{i})))-2;
    for j=1:nstacks
        nimages = length(dir(sprintf('%s/%s/%s_%d', image_dir, folders{i}, folders{i}, j)))-2;
        image_stack_tmp = zeros(npixels, npixels, nchannels, nimages, 'uint8');
        total_intensity = zeros(nimages, 1);
        for k=1:nimages
            im_tmp = imread(sprintf('%s/%s/%s_%d/%s_%d-%d.tif', image_dir, folders{i}, folders{i}, j, folders{i}, j, 10000+k-1));
            image_stack_tmp(:,:,:,k) = imresize(immultiply(im_tmp, 1.5), [npixels npixels]);
            total_intensity(k) = sum(double(im_tmp(:)));
        end
%         [~, idx] = max(squeeze(sum(sum(sum(double(image_stack_tmp), 1), 2), 3)));
        [~, idx] = max(total_intensity);
        
        image_stacks = cat(5, image_stacks, image_stack_tmp(:,:,:,idx-nimages_stack_UD:idx+nimages_stack_UD));
        time = [time time_for_folders(i)];
    end
end

nimages = length(time);

%%

for i=[1:7 10 15 17 18 22 25:30 32 33 35 41:46]
    image_stacks(:,:,:,:,i) = image_stacks(end:-1:1,:,:,:,i);
end

for i=[1 2 7 9 11 13 16 22 25 29:31 33 36 37 39 40 42 44]
    image_stacks(:,:,:,:,i) = image_stacks(:,end:-1:1,:,:,i);
end

%%

% number of angles used to compute the pairwise alignments
nrot = 36;
% do not consider shifts
nshifts = 0;
shift_max = 0;
% show progress bar
show_waitbar = true;

images_tmp = uint8(squeeze(mean(double(image_stacks(:,:,:,:,:)), 4)));
% images_tmp = squeeze(max(image_stacks(:,:,:,:,:), [],4));
% images_tmp(:,:,1,:,:) = immultiply(images_tmp(:,:,1,:,:), 0.5);
[R, W] = compute_pairwise_alignments(images_tmp, nrot, shift_max, nshifts, show_waitbar);

%%
% compute optimal rotations + embedding coordinates using vector diffusion maps
% use default kernel scale
eps = median(W(:));
% only compute first embedding coordinate
ncomps = 1;

[R_opt, embed_coord, D2] = vdm(R, W, eps, ncomps);

% register images using optimal rotations
image_stacks2 = register_all_images(image_stacks, R_opt);

%%
figure;
for i=1:nimages
    subplot(6,8,i)
    imshow(images_tmp(:,:,:,i));
end

figure;
for i=1:nimages
    subplot(6,8,i)
    imshow((image_stacks(:,:,:,nimages_stack_UD+1,i)));
end

figure;
for i=1:nimages
    subplot(6,8,i)
    imshow((image_stacks2(:,:,:,nimages_stack_UD+1,i)));
end


%%
W = zeros(nimages);

% dist_tmp = zeros(4, 1);
% for i=1:nimages
%     im1 = image_stacks(:,:,:,:,i);
%     for j=1:i-1
%         im2 = image_stacks(:,:,:,:,j);
%         dist_tmp(1) = sum((im1(:)-im2(:)).^2);
%
%         im2 = image_stacks(end:-1:1,:,:,:,j);
%         dist_tmp(2) = sum((im1(:)-im2(:)).^2);
%
%         im2 = image_stacks(:,end:-1:1,:,:,j);
%         dist_tmp(3) = sum((im1(:)-im2(:)).^2);
%
%         im2 = image_stacks(end:-1:1,end:-1:1,:,:,j);
%         dist_tmp(4) = sum((im1(:)-im2(:)).^2);
%
%         W(i,j) = min(dist_tmp);
%     end
% end

for i=1:nimages
    im1 = image_stacks2(:,:,:,:,i);
    for j=1:i-1
        im2 = image_stacks2(:,:,:,:,j);
        
        W(i,j) = sum((im1(:)-im2(:)).^2);
    end
end

W = W + W';

%%

eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

figure;
plot(V(:,2), time, '.')

[~, I] = sort(V(:,2));
figure;
for i=1:nimages
    subplot(6,8,i)
    imshow((image_stacks2(:,:,:,nimages_stack_UD+1,I(i))));
end


