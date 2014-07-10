clear all
close all

% %% read in movies
% 
% npixels = 100;
% channel_his = 1;
% channel_sqh = 2;
% 
% file_names = {'bomyi_emb01_gast01.avi'; 'bomyi_emb02_gast02.avi';...
%     '14_0623/emb01_cell_gast.avi'; '14_0623/emb02_cell_gast.avi'; ...
%     '14_0624/emb01_cell_gast.avi'; '14_0624/emb02_cell_gast.avi'};
% 
% nmovies = length(file_names);
% images = [];
% time = [];
% movie_idx = [];
% 
% % movie_start_idx = [108 62 100 101];
% % movie_start_idx = [103 66 92 93];
% 
% % using when dent forms in sqh
% movie_start_idx = [130 86 108 115];
% 
% figure;
% for i=1:nmovies
%     
%     [images_tmp, times_tmp] = read_video(file_names{i}, npixels);
%     images_tmp = images_tmp(:, :, channel_his, :);
%     images_tmp = squeeze(images_tmp);
%     
%     total_intensity_sqh = sum(sum(images_tmp));
% %     total_intensity_sqh = max(max(images_tmp));
%     
% %     [~, idx] = max(squeeze(total_intensity_sqh));
% %     idx = movie_start_idx(i)-20;
%     
%     [~, idx] = min(squeeze(total_intensity_sqh));
%     idx = idx - 25;
%     
%     subplot(3,2,i);
%     plot(times_tmp, squeeze(total_intensity_sqh),'.')
%     hold on
%     plot(times_tmp(idx), total_intensity_sqh(idx),'.r')
%     
%     
%     [images_tmp, times_tmp] = read_video(file_names{i}, npixels);
%     images_tmp = images_tmp(:, :, channel_his, :);
%     images_tmp = squeeze(images_tmp);
%     
%     times_tmp = times_tmp(idx:end) - times_tmp(idx);
%     images_tmp = images_tmp(:, :, idx:end);
%     
%     images = cat(3, images, images_tmp);
%     time = [time; times_tmp];
%     movie_idx = [movie_idx; i*ones(size(times_tmp))];
%     
% end
% 
% 
% %% adjust images
% 
% % H = fspecial('disk',5);
% % image_fn = @(im1) imfilter(imadjust(im1),H,'replicate');
% 
% % image_fn = @(im1) imadjust(im1);
% 
% image_fn = @(im1) uint8(double(im1) - mean(double(im1(:))));
% 
% %% compute scattering coefficients
% 
% addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
% % addpath '../../../MATLAB/scatnet-0.2';
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
% x = double(image_fn(images(:,:,1)));
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
%     x = double(image_fn(images(:,:,i)));
%     
%     Sx = scat(x, Wop);
%     
%     sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
% end
% 
% save('movie_data_zero.mat','images','time','movie_idx','sx_all','file_names','nmovies');

% return

%% load data

load('movie_data_zero.mat');

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

figure;
bar(diag(D(1:10,1:10))/sum(diag(D)))

figure;
semilogy(diag(D(1:10,1:10))/sum(diag(D)))

train_coeff = train_sx_all * V;
test_coeff = test_sx_all * V;
% all_coeff = (sx_all - repmat(mean_data, length(time), 1)) * V;


%%

for i=1:2
    figure;
    scatter(train_times,train_coeff(:,i), 200, train_movie_idx,'.')
    hold on
    scatter(test_times,test_coeff(:,i), 200, test_movie_idx,'.')
end

%% plot PCA projection

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times, '.')

figure;
scatter3(train_coeff(:,1), train_coeff(:, 2), train_coeff(:, 3),50, train_times, '.')

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
scatter(test_coeff(:,1), test_coeff(:, 2), 50, 'k','.')

%% kernel interpolation

PCA_idx = 1:4;

% weights = pdist2(test_coeff(:, PCA_idx), train_coeff(:, PCA_idx));
% kernel_eps = median(weights(:))/5;
% weights = exp(-weights.^2 / kernel_eps^2);
% for i=1:length(test_times)
%     weights(i,:) = weights(i,:) / sum(weights(i,:));
% end
%
% interp_test_times = weights * train_times;
%
% figure;
% plot(test_times, interp_test_times,'.')
% hold on
% plot([min(time) max(time)],[min(time) max(time)],'-r')
% xlabel('true times')
% ylabel('predicted times')



%% kernel ridge regression

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

%%

h1 = figure;
h2 = figure;
h3 = figure;

for j=1:nmovies
    test_movies = j;
    train_movies = setdiff(1:nmovies, test_movies);
    
    idx = ismember(movie_idx, train_movies);
    train_sx_all = sx_all(idx, :);
    train_times = time(idx);
    train_movie_idx = movie_idx(idx);
    
    idx = ismember(movie_idx, test_movies);
    test_sx_all = sx_all(idx, :);
    test_times = time(idx);
    test_movie_idx = movie_idx(idx);
    
    % mean-center
    mean_data = mean(train_sx_all);
    
    train_sx_all = train_sx_all - repmat(mean_data, length(train_times), 1);
    test_sx_all = test_sx_all - repmat(mean_data, length(test_times), 1);
    
    % PCA
    [V, D] = PCA(train_sx_all, 100);

    train_coeff = train_sx_all * V;
    test_coeff = test_sx_all * V;
    

    % kernel ridge regression
    
    lambda = logspace(-4, 2, 100);
    err = zeros(size(lambda));
    
    K = squareform(pdist(train_coeff(:, PCA_idx))).^2;
    k = pdist2(train_coeff(:, PCA_idx), test_coeff(:, PCA_idx)).^2;
    eps = median(K(:));
    
    K = exp(-K/eps);
    k = exp(-k/eps);
    
    for i=1:length(lambda)
        
        interp_test_times = train_times' * inv(K + lambda(i) * eye(size(K,1))) * k;
        
        err(i) = sum((test_times - interp_test_times').^2);
    end
    
    [~, lambda_idx] = min(err);
    interp_test_times = train_times' * inv(K + lambda(lambda_idx) * eye(size(K,1))) * k;
    
    figure(h1);
    subplot(3,2,j);
    plot(test_times, interp_test_times', '.')
    hold on
    plot([min(test_times) max(test_times)],[min(test_times) max(test_times)],'-r')
    xlabel('true times')
    ylabel('predicted times')
        title(sprintf('mean abs error = %2.1f', mean(abs(test_times-interp_test_times'))))

    figure(h2);
    subplot(3,2,j);
    plot(test_times, abs(test_times-interp_test_times'), '.')
    xlabel('true times')
    ylabel('abs error')
    
    figure(h3);
    subplot(3,2,j);
    plot(test_times, interp_test_times', '.')
    hold on
    idx = find(abs(test_times-interp_test_times') > 10);
        plot(test_times(idx), interp_test_times(idx)', '.r')

    plot([min(test_times) max(test_times)],[min(test_times) max(test_times)],'-r')
    xlabel('true times')
    ylabel('predicted times')
    
    mean(abs(test_times-interp_test_times'))
end







