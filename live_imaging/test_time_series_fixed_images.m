clear all
close all


%% parameters

% number of pixels to subsample the movies to (doing the scattering
% transform on the original high-resolution movies takes too long)
npixels = 100;

% channel of the movie which contains the revevant signal (red = 1, green =
% 2, blue = 3)
channel = 1;

% path where scattering transform code is stored
% scatnet_path = 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
scatnet_path = '../../../MATLAB/scatnet-0.2';

% directory name where images are stored
dir_name = '../membrane_pictures/14_0501_dpERK_late';

nimages = 132;

%% define function to adjust each image
image_fn = @(im1) imadjust(im1);


%% read in images

nchannels = 3;
fixed_images = zeros(npixels, npixels, nchannels, nimages, 'uint8');

for i=1:nimages
    A = imread(sprintf('%s/emb%02d.tif', dir_name, i));
    A = imresize(A, [npixels npixels]);
    A(:, :, channel) = image_fn(A(:,:,channel));
    
    fixed_images(:, :, :, i) = A;
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
x = double(fixed_images(:,:,1,1));
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
sx_all = zeros(nimages, size(Sx_mat, 1), nchannels);

% compute scattering invariants for each image
for i=1:nimages
    i
    for j=1:nchannels
        x = double(fixed_images(:,:,j,i));

        Sx = scat(x, Wop);

        sx_all(i, :, j) = mean(mean(format_scat(Sx),2),3)';
    end
end

return


%% read in movies

% number of movies
nmovies = length(file_names);

% array which stores images
images = [];

% array which stores time for each image
time = [];

% array which stores the index of the movie which generated the
% corresponding image frame
movie_idx = [];

for i=1:nmovies

    [images_tmp, times_tmp] = read_video(sprintf('%s',file_names{i}), npixels);
    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);
    
    images = cat(3, images, images_tmp);
    

    time = [time; times_tmp];
    movie_idx = [movie_idx; i*ones(size(times_tmp))];

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
sx_all = zeros(length(time), size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:length(time)
    i
    x = double(image_fn(images(:,:,i)));

    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%% mean-center data (for PCA)

mean_data = mean(sx_all);

sx_all = sx_all - repmat(mean_data, length(time), 1);

%% PCA
[V, D] = PCA(sx_all, 100);

figure;
bar(diag(D(1:10,1:10))/sum(diag(D)))
xlabel('k')
ylabel('\lambda_k')
title('eigenvalues from PCA of scattering transform coefficients')

%% project onto PCA basis

all_coeff = sx_all * V;

%% plot PCA projection

figure;
scatter(all_coeff(:,1), all_coeff(:, 2), 50, time, '.')
colorbar
xlabel('projection onto PC 1')
ylabel('projection onto PC 2')
title('Projection of scattering coefficients onto first two PCs; data colored by time')

%% kernel ridge regression for time prediction

% projections to use in regression
PCA_idx = 1:5;

idx1 = find(movie_idx == 1);
idx2 = find(movie_idx == 2);

W = pdist2(all_coeff(idx1, PCA_idx), all_coeff(idx2, PCA_idx));

figure;
imagesc(-W)

figure; 
imagesc(-log(W))