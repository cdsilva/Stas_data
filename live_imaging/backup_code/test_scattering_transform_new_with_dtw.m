clear all
close all

%%

image_fn = @(image) adapthisteq(crop_image(image, 200, 0));

scatnet_path = 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
% scatnet_path = '../../../MATLAB/scatnet-0.2';

channel = 1;
npixels = 100;

[images, time] = read_video('14_0624/emb01_hisRFP_gastrulation.avi', npixels);
images = images(:, :, channel, :);
images = squeeze(images);

% load('../image_analysis_paper/late_images_aligned.mat','image_set_raw_ordered');
% time = (1:size(image_set_raw_ordered, 4))';
% images = zeros(npixels, npixels, size(image_set_raw_ordered, 4), 'uint8');
% for i=1:length(time)
%     images(:,:,i) = imresize(image_set_raw_ordered(:,:,channel,i), [npixels npixels]);
% end

nimages = length(time);
subplot_dim1 = ceil(sqrt(nimages));
subplot_dim2 = ceil(nimages / subplot_dim1);

%%

[images_movie, times_movie] = read_video('bomyi_emb01_gast01.avi', npixels);
images_movie = images_movie(:, :, channel, :);
images_movie = squeeze(images_movie);

nimages_movie = length(times_movie);
subplot_dim1_movie = ceil(sqrt(nimages_movie));
subplot_dim2_movie = ceil(nimages_movie / subplot_dim1_movie);

%%

figure;
for i=1:nimages
    images(:,:,i) = image_fn(images(:,:,i));
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(images(:,:,i));
end

figure;
for i=1:nimages_movie
    images_movie(:,:,i) = image_fn(images_movie(:,:,i));
    
    subplot(subplot_dim1_movie, subplot_dim2_movie, i)
    imshow(images_movie(:,:,i));
end


%% compute scattering coefficients

addpath(scatnet_path);
addpath_scatnet

sx_tmp = compute_scattering_coeff(images(:,:,1));

sx_all = zeros(length(nimages), length(sx_tmp));
sx_all_movie = zeros(length(nimages_movie), length(sx_tmp));

% compute scattering invariants for each image
H = waitbar(0, 'Computing scattering coefficients');
for i=1:nimages
	waitbar(i/(nimages+nimages_movie), H);
    sx_all(i,:) = compute_scattering_coeff(images(:,:,i));
end

for i=1:nimages_movie
	waitbar((i+nimages)/(nimages+nimages_movie), H);
    sx_all_movie(i,:) = compute_scattering_coeff(images_movie(:,:,i));
end
close(H);

%%


[V, D] = PCA(sx_all, 10);

figure;
bar(diag(D) / sum(diag(D)))

PCA_idx = 1:3;

%%

% W = pdist2(sx_all(randperm(nimages), PCA_idx), sx_all_movie(:, PCA_idx));
W = pdist2(sx_all(:, PCA_idx), sx_all_movie(:, PCA_idx));

figure;
imagesc(W)


%%

[D_dtw, warp_path] = DTW(W);

figure;
imagesc(D_dtw)

figure;
imagesc(warp_path)
colormap(gray)

% sum(sum(D_dtw.*warp_path))
D_dtw(end, end)
