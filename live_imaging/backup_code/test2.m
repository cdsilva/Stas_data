clear all
close all

%% read in movie

train_file_name = 'bomyi_emb01_gast01.avi';

test_file_name = 'bomyi_emb02_gast02.avi';

npixels = 100;
channel = 1;

[train_images, train_times] = read_video(train_file_name, npixels);
train_images = train_images(:, :, channel, :);


[test_images, test_times] = read_video(test_file_name, npixels);
test_images = test_images(:, :, channel, :);


%%
% addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
addpath '../../../MATLAB/scatnet-0.2';
addpath_scatnet

%% compute scattering coefficients
% set scattering transform parameters
filt_opt = struct();
filt_rot_opt = struct();
% oversampling is set to infty
% so that scattering of original and rotated
% image will be sampled at exactly the same points
scat_opt = struct();
scat_opt.oversampling = 10;

% compute scattering transform of first image
x = double(train_images(:,:,1));
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
train_sx_all = zeros(length(train_times), size(Sx_mat, 1));
test_sx_all = zeros(length(test_times), size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:length(train_times)
    i
    x = double(train_images(:,:,i));
    
    Sx = scat(x, Wop);
    
    train_sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

for i=1:length(test_times)
    i
    x = double(test_images(:,:,i));
    
    Sx = scat(x, Wop);
    
    test_sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%% mean-center

train_sx_all = train_sx_all - repmat(mean(train_sx_all), length(train_times), 1);
test_sx_all = test_sx_all - repmat(mean(test_sx_all), length(test_times), 1);

%% PCA
[V, D] = PCA(train_sx_all, 2);

train_coeff = train_sx_all * V;
test_coeff = test_sx_all * V;

%% plot PCA projection
figure;
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition', [0 0 8 8])
scatter(train_coeff(:,1), train_coeff(:, 2), 200, train_times, '.')
hold on
ind = 35:57;
scatter(train_coeff(ind,1), train_coeff(ind, 2), 200, train_times(ind), 'o')
plot(test_coeff(:,1),test_coeff(:,2),'xk')
xlabel('scattering transform coefficients projected onto first PC')
ylabel('scattering transform coefficients projected onto second PC')
legend('training data', 'testing data')

%% kernel interpolation

kernel_dist = pdist2(test_coeff, train_coeff).^2;
eps = median(kernel_dist(:))/10;
kernel_dist = exp(-kernel_dist/eps);

for i=1:length(test_times)
    kernel_dist(i,:) = kernel_dist(i,:) / sum(kernel_dist(i,:));
end

interp_test_times = kernel_dist * train_times;

figure;
plot(test_times, interp_test_times, '.')

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
scatter(test_coeff(:,1),test_coeff(:,2),50, interp_test_times, '.')

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
scatter(test_coeff(:,1),test_coeff(:,2),50, test_times, '.')

%% fit linear model

B = regress(train_times,[ones(length(train_times), 1) train_coeff train_coeff.^2]);

figure;
plot([ones(length(train_times), 1) train_coeff train_coeff.^2] * B, train_times, '.')

figure;
plot([ones(length(test_times), 1) test_coeff test_coeff.^2] * B, test_times, '.')

%% DMAPS

W = squareform(pdist(train_coeff)).^2;
eps = median(W(:))/10;
[V_dmaps, D_dmaps] = dmaps(W, eps, 10);

figure;
plot(V_dmaps(:,2),train_times,'.')

corr(V_dmaps(:,2),train_times, 'type', 'spearman')

figure;
scatter(V_dmaps(:,2),V_dmaps(:,4),50, train_times)

%%
W = squareform(pdist([train_coeff; test_coeff])).^2;
eps = median(W(:))/30;
[V_dmaps, D_dmaps] = dmaps(W, eps, 10);

figure;
plot(V_dmaps(:,2),[train_times; test_times],'.')


%%

Wuu = squareform(pdist(test_coeff)).^2;
Wul = pdist2(test_coeff, train_coeff).^2;

eps = median(Wuu(:));
Wuu = exp(-Wuu/eps);
Wul = exp(-Wul/eps);

fl = train_times - mean(train_times);
Du = diag(sum(Wuu));

fu = inv(Du - Wuu) * Wul * fl;

figure;
plot(test_times, fu, '.')




