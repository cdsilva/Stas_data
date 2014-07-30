clear all
close all


%% parameters

% number of pixels to subsample the movies to (doing the scattering
% transform on the original high-resolution movies takes too long)
npixels = 100;

% channel of the movie which contains the revevant signal (red = 1, green =
% 2, blue = 3)
channel = 1;

% directory where movies are stored (relative to the current directory)
file_dir = 'grayscale_processed_movies';

% file names of movies
file_names = {'0623_emb01_gray.avi'; '0623_emb02_gray.avi'; ...
    '0709_emb02_gray.avi'; 'bomyi_emb01_gray.avi'; 'bomyi_emb02_gray.avi';
    'emb01_cell_gast_gray.avi'; 'emb02_cell_gast_gray.avi'};

% start frame for each movie (when the "dent" begins)
frame_start = [126; 71; 107; 1; 1; 93; 99];

% which movies to use for training
train_movies = 1:6;
% which movies to use for testing
test_movies = [7];

% path where scattering transform code is stored
scatnet_path = 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';

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

    [images_tmp, times_tmp] = read_video(sprintf('%s/%s', file_dir,file_names{i}), npixels);
    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);
    images_tmp = images_tmp(:, :, frame_start(i):end);
    
    images = cat(3, images, images_tmp);
    
    times_tmp = times_tmp(frame_start(i):end) - times_tmp(frame_start(i));

    time = [time; times_tmp];
    movie_idx = [movie_idx; i*ones(size(times_tmp))];

end

%% define function to adjust each image
image_fn = @(im1) imadjust(im1);

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

%% select training and testing data

idx = ismember(movie_idx, train_movies);
% scattering coefficients for training data
train_sx_all = sx_all(idx, :);
% times for training data
train_times = time(idx);
% movie idx for training data
train_movie_idx = movie_idx(idx);

idx = ismember(movie_idx, test_movies);
% scattering coefficients for test data
test_sx_all = sx_all(idx, :);
% times for test data
test_times = time(idx);
% movie idx for test data
test_movie_idx = movie_idx(idx);

%% mean-center data (for PCA)

mean_data = mean(train_sx_all);

train_sx_all = train_sx_all - repmat(mean_data, length(train_times), 1);
test_sx_all = test_sx_all - repmat(mean_data, length(test_times), 1);

%% PCA
[V, D] = PCA(train_sx_all, 100);

figure;
bar(diag(D(1:10,1:10))/sum(diag(D)))
xlabel('k')
ylabel('\lambda_k')
title('eigenvalues from PCA of scattering transform coefficients')

%% project onto PCA basis

train_coeff = train_sx_all * V;
test_coeff = test_sx_all * V;
all_coeff = (sx_all - repmat(mean_data, length(time), 1)) * V;

%% plot PCA projection

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times, '.')
colorbar
xlabel('projection onto PC 1')
ylabel('projection onto PC 2')
title('Projection of scattering coefficients onto first two PCs; data colored by time')

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
colorbar
hold on
scatter(test_coeff(:,1), test_coeff(:, 2), 50, 'k','.')
legend('training data','testing data','location','best')
xlabel('projection onto PC 1')
ylabel('projection onto PC 2')
title('Projection of scattering coefficients onto first two PCs; data colored by time')

%% kernel ridge regression for time prediction

% projections to use in regression
PCA_idx = 1:4;

% find optimal value of regularization parameter lambda by iterating over
% many values and computing error
lambda = logspace(-4, 2, 100);
err = zeros(size(lambda));

% distance matrices
K = squareform(pdist(train_coeff(:, PCA_idx))).^2;
k = pdist2(train_coeff(:, PCA_idx), test_coeff(:, PCA_idx)).^2;

% kernel scale
eps = median(K(:));

% compute kernel functin of distances
K = exp(-K/eps);
k = exp(-k/eps);

% calculate error in predictions for each value of lambda
for i=1:length(lambda)
    
    interp_test_times = train_times' * inv(K + lambda(i) * eye(size(K,1))) * k;
    
    err(i) = sum((test_times - interp_test_times').^2);
end

% figure;
% semilogx(lambda, err, '.')
% xlabel('\lambda')
% ylabel('error in predicted times of test data')

% use value of lambda that yields minimum error 
[~, lambda_idx] = min(err);
interp_test_times = train_times' * inv(K + lambda(lambda_idx) * eye(size(K,1))) * k;

% plot predicted times for testing data
figure;
scatter(test_times, interp_test_times',50, test_movie_idx, '.')
colormap(cool)
hold on
plot([min(test_times) max(test_times)],[min(test_times) max(test_times)],'-r')
xlabel('true times')
ylabel('predicted times using kernel ridge regression')
