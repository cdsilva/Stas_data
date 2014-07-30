clear all
close all

%%

image_fn = @(image) adapthisteq(medfilt2(imresize(image, [100 100]), [5 5]),'distribution','exponential');

scatnet_path = 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';

channel = 1;
npixels = 100;

[images, time] = read_video('14_0624/emb02_hisRFP_gastrulation.avi', npixels);
images = images(:, :, channel, :);
images = squeeze(images);

nimages = length(time);
subplot_dim1 = ceil(sqrt(nimages));
subplot_dim2 = ceil(nimages / subplot_dim1);

figure;
for i=1:nimages
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(image_fn(images(:,:,i)));
end

%%

[images_movie, times_movie] = read_video('bomyi_emb01_gast01.avi', npixels);
images_movie = images_movie(:, :, channel, :);
images_movie = squeeze(images_movie);

nimages_movie = length(times_movie);
subplot_dim1_movie = ceil(sqrt(nimages_movie));
subplot_dim2_movie = ceil(nimages_movie / subplot_dim1_movie);

figure;
for i=1:nimages_movie
    subplot(subplot_dim1_movie, subplot_dim2_movie, i)
    imshow(image_fn(images_movie(:,:,i)));
end

%%

figure;
for i=1:nimages
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(image_fn(images(:,:,i)));
end

figure;
for i=1:nimages_movie
    subplot(subplot_dim1_movie, subplot_dim2_movie, i)
    imshow(image_fn(images_movie(:,:,i)));
end


%% compute scattering coefficients

addpath(scatnet_path);
addpath_scatnet

% set scattering transform parameters
filt_opt = struct();
filt_rot_opt = struct();
% oversampling is set to infty
% so that scattering of original and rotated
% image will be sampled at exactly the same points
scat_opt = struct();
scat_opt.oversampling = 10;

% compute scattering transform of first image
x = double(image_fn(images(:,:,1)));
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
sx_all = zeros(length(nimages), size(Sx_mat, 1));
sx_all_movie = zeros(length(nimages_movie), size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:nimages
    i
    x = double(image_fn(images(:,:,i)));

    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

for i=1:nimages_movie
    i
    x = double(image_fn(images_movie(:,:,i)));

    Sx = scat(x, Wop);

    sx_all_movie(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%%


[V, D] = PCA(sx_all, 10);

figure; 
scatter(sx_all_movie*V(:,1),sx_all_movie*V(:,2),50, 1:nimages_movie,'.')

figure; 
subplot(1,2,1)
scatter((sx_all-repmat(mean(sx_all), nimages, 1))*V(:,1),(sx_all-repmat(mean(sx_all), nimages, 1))*V(:,2),50, 1:nimages,'o')
subplot(1,2,2)
scatter((sx_all_movie-repmat(mean(sx_all_movie), nimages_movie, 1))*V(:,1),(sx_all_movie-repmat(mean(sx_all_movie), nimages_movie, 1))*V(:,2),50, 1:nimages_movie,'.')

figure; 
scatter((sx_all-repmat(mean(sx_all), nimages, 1))*V(:,1),(sx_all-repmat(mean(sx_all), nimages, 1))*V(:,2),50, 1:nimages,'o')
hold on
scatter((sx_all_movie-repmat(mean(sx_all_movie), nimages_movie, 1))*V(:,1),(sx_all_movie-repmat(mean(sx_all_movie), nimages_movie, 1))*V(:,2),50, 1:nimages_movie,'.')

%%

W = pdist2(sx_all, sx_all_movie);

figure; 
imagesc(W)

figure; 
imagesc(-log(W))