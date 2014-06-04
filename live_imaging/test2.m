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

[V, D] = PCA(train_sx_all, 2);

train_coeff = train_sx_all * V;
test_coeff = test_sx_all * V;

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
plot(test_coeff(:,1),test_coeff(:,2),'.')

%%

B = regress(train_times,[ones(length(train_times), 1) train_coeff]);

figure;
plot([ones(length(train_times), 1) train_coeff] * B, train_times, '.')

coeff2 = test_sx_all * V;
figure;
plot([ones(length(test_times), 1) coeff2] * B, test_times, '.')

%%

W = squareform(pdist(train_coeff)).^2;
eps = median(W(:));
[V_dmaps, D_dmaps] = dmaps(W, eps/10, 10);

figure;
plot(V_dmaps(:,2),train_times,'.')

corr(V_dmaps(:,2),train_times, 'type', 'spearman')

figure;
scatter(V_dmaps(:,2),V_dmaps(:,4),50, train_times)



