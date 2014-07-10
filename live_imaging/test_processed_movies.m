clear all
close all

% %% read in movies
% 
% npixels = 100;
% channel = 1;
% 
% file_dir = 'grayscale_processed_movies';
% 
% file_names = {'0623_emb01_gray.avi'; '0623_emb02_gray.avi'; ...
%     '0709_emb02_gray.avi'; 'bomyi_emb01_gray.avi'; 'bomyi_emb02_gray.avi';
%     'emb01_cell_gast_gray.avi'; 'emb02_cell_gast_gray.avi'};
% 
% frame_start = [126; 71; 107; 1; 1; 93; 99];
% 
% nmovies = length(file_names);
% images = [];
% time = [];
% movie_idx = [];
% 
% for i=1:nmovies
% 
%     [images_tmp, times_tmp] = read_video(sprintf('%s/%s', file_dir,file_names{i}), npixels);
%     images_tmp = images_tmp(:, :, channel, :);
%     images_tmp = squeeze(images_tmp);
%     images_tmp = images_tmp(:, :, frame_start(i):end);
%     
%     images = cat(3, images, images_tmp);
%     
%     times_tmp = times_tmp(frame_start(i):end) - times_tmp(frame_start(i));
%     
%     time = [time; times_tmp];
%     movie_idx = [movie_idx; i*ones(size(times_tmp))];
% 
% end
% 
% %% compute scattering coefficients
% 
% % addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
% addpath '../../../MATLAB/scatnet-0.2';
% addpath_scatnet
% 
% % set scattering transform parameters
% filt_opt = struct();
% filt_rot_opt = struct();
% % oversampling is set to infty
% % so that scattering of original and rotated
% % image will be sampled at exactly the same points
% scat_opt = struct();
% scat_opt.oversampling = 10;
% 
% % compute scattering transform of first image
% x = double(images(:,:,1));
% % define wavelet transforms
% Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
% Sx = scat(x, Wop);
% Sx_mat = format_scat(Sx);
% 
% % store scattering invariants (one row per image)
% sx_all = zeros(length(time), size(Sx_mat, 1));
% 
% % compute scattering invariants for each image
% for i=1:length(time)
%     i
%     x = double(images(:,:,i));
% 
%     Sx = scat(x, Wop);
% 
%     sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
% end
% 
% save('movie_data_grayscale.mat','images','time','movie_idx','sx_all','file_names','nmovies');
% 
% return

%% load data

load('movie_data_grayscale.mat');


%% select data
train_movies = [1,3,5];
test_movies = setdiff(1:nmovies, train_movies);

idx = ismember(movie_idx, train_movies);
train_sx_all = sx_all(idx, :);
train_times = time(idx);
train_movie_idx = movie_idx(idx);

idx = ismember(movie_idx, test_movies);
test_sx_all = sx_all(idx, :);
test_times = time(idx);
test_movie_idx = movie_idx(idx);

%% mean-center
mean_data = mean(train_sx_all);

train_sx_all = train_sx_all - repmat(mean_data, length(train_times), 1);
test_sx_all = test_sx_all - repmat(mean_data, length(test_times), 1);

%% PCA
[V, D] = PCA(train_sx_all, 100);

if mean(V(:,1)) > 0
    V(:,1) = -V(:,1);
end

figure;
bar(diag(D(1:10,1:10))/sum(diag(D)))

train_coeff = train_sx_all * V;
test_coeff = test_sx_all * V;
all_coeff = (sx_all - repmat(mean_data, length(time), 1)) * V;

%% plot PCA projection

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times, '.')

figure;
scatter3(train_coeff(:,1), train_coeff(:, 2), train_coeff(:, 3),50, train_times, '.')

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
scatter(test_coeff(:,1), test_coeff(:, 2), 50, 'k','.')



%% kernel ridge regression

PCA_idx = 1:4;

lambda = logspace(-4, 2, 100);
err = zeros(size(lambda));

K = squareform(pdist(train_coeff(:, PCA_idx))).^2;
k = pdist2(train_coeff(:, PCA_idx), test_coeff(:, PCA_idx)).^2;
eps = median(K(:));

K = exp(-K/eps);
k = exp(-k/eps);

for i=1:length(lambda)
    %     lambda = 0.1;
    
    interp_test_times = train_times' * inv(K + lambda(i) * eye(size(K,1))) * k;
    
    % figure;
    % plot(test_times, interp_test_times','.')
    % hold on
    % plot([0 max(time)],[0 max(time)],'-r')
    % xlabel('true times')
    % ylabel('predicted times')
    
    err(i) = sum((test_times - interp_test_times').^2);
end

figure;
semilogx(lambda, err, '.')

[~, lambda_idx] = min(err);
interp_test_times = train_times' * inv(K + lambda(lambda_idx) * eye(size(K,1))) * k;

figure;
scatter(test_times, interp_test_times',50, test_movie_idx, '.')
colormap(cool)
hold on
plot([min(test_times) max(test_times)],[min(test_times) max(test_times)],'-r')
xlabel('true times')
ylabel('predicted times')







