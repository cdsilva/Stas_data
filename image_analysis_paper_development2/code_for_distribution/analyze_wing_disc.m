clear all
close all

%% read in images

% number of pixels
npixels = 100;

% image directory
image_dir = 'wing_disc';

image_files = dir('wing_disc');

% number of images
nimages = length(image_files)-2;

% image prefix
% image_name = 'emb';
% channel where nuclear signal is stored (1=red, 2=green, 3=blue)
% nuclear_channel = 1;
% number of channels in images
nchannels = 3;

% process images
images = zeros(npixels, npixels, nchannels, nimages, 'uint8');
time = zeros(nimages, 1);
for i=1:nimages
    
    % read image
    im_tmp = imread(sprintf('%s/%s', image_dir, image_files(i+2).name));
    
    % extract relevant channel; resize to number of pixels
    im_tmp = imresize(im_tmp, [npixels npixels]);
    
    % crop image
    %     im_tmp = crop_image(im_tmp, nuclear_channel);
    
    % normalize+scale nuclear signal
    %     im_tmp(:,:,nuclear_channel) = normalize_nuclear_signal(im_tmp(:,:,nuclear_channel));
    %     im_tmp(:, :, nuclear_channel) = immultiply(im_tmp(:, :, nuclear_channel), 0.5);
    
    im_tmp = immultiply(im_tmp, 1.5);
    
    images(:, :, :, i) = im_tmp;
    
    k = strfind(image_files(i+2).name, '_');
    time_str = image_files(i+2).name(k(end-1)+1:k(end)-1);
    time_str = strrep(time_str, ',', '.');
    k = strfind(time_str, '-');
    time1 = str2double(time_str(1:k-1));
    time2 = str2double(time_str(k+1:end));
    
    time(i) = (time1 + time2)/2;
    
    sample_idx = str2double(image_files(i+2).name(end-4:end-4));
    
    if time1 == 72
        if any(sample_idx == [1 2 7 8])
            images(:, :, :, i) = images(end:-1:1, :, :, i);
        end
        
    elseif time1 == 76.5
        if any(sample_idx == [3 5])
            images(:, :, :, i) = images(end:-1:1, :, :, i);
        end
    elseif time1 == 79
        if any(sample_idx == [3:8])
            images(:, :, :, i) = images(end:-1:1, :, :, i);
        end
    elseif time1 == 89
        if any(sample_idx == [1 5:7 9])
            images(:, :, :, i) = images(end:-1:1, :, :, i);
        end
    elseif time1 == 100
        if any(sample_idx == [1 5])
            images(:, :, :, i) = images(end:-1:1, :, :, i);
        end
    elseif time1 == 110.5
        if any(sample_idx == [2 4])
            images(:, :, :, i) = images(end:-1:1, :, :, i);
        end
    else
        disp('INVALID TIME')
        return
    end
    
    
end


%%

% idx = setdiff(1:nimages, [5 10 1 7 4 2 6 13 8 11 12 9]);
% 
% time = time(idx);
% images = images(:,:,:,idx);
% nimages = length(idx);

%% show raw images

figure;
for i=1:nimages
    subplot(6, 8, i)
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
eps = median(W(:));
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
    subplot(6, 8, i)
    imshow(images_registered_ordered(:,:,:,i));
end

%% compare ranks
figure;
plot(compute_ranks(time),compute_ranks(embed_coord(:,1)),'.')
xlabel('expert rank')
ylabel('algorithm rank')

%%
figure;
plot(time,compute_ranks(embed_coord(:,1)),'.')
xlabel('time')
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
