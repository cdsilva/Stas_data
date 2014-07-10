clear all
close all

%% read in movie

train_file_name = 'bomyi_emb01_gast01.avi';

test_file_name = 'bomyi_emb02_gast02.avi';

npixels = 100;
channel = 1;

[train_images, train_times] = read_video(train_file_name, npixels);
train_images = train_images(:, :, channel, :);
train_images = squeeze(train_images);

[test_images, test_times] = read_video(test_file_name, npixels);
test_images = test_images(:, :, channel, :);
test_images = squeeze(test_images);

test_dir_names = {'../membrane_pictures/14_0403_dpERK_late'; '../membrane_pictures/14_0501_dpERK_late'};
num_images = [46; 132];

for j=1:length(num_images)
    for i=1:num_images(j)
        A = imread(sprintf('%s/emb%02d.tif', test_dir_names{j}, i));
        A = A(:, :, channel);
        A = imresize(A, [npixels, npixels]);
        test_images = cat(3, test_images, A);
        test_times = [test_times; 0];
    end
end

%% adjust images

H = fspecial('disk',5);

for i=1:length(train_times)
    train_images(:, :, i) = imadjust(train_images(:, :, i));
    train_images(:, :, i) = imfilter(train_images(:, :, i),H,'replicate');
end

for i=1:length(test_times)
    test_images(:, :, i) = imadjust(test_images(:, :, i));
    test_images(:, :, i) = imfilter(test_images(:, :, i),H,'replicate');
    
end

%%
addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
% addpath '../../../MATLAB/scatnet-0.2';
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

mean_data = mean(train_sx_all);

train_sx_all = train_sx_all - repmat(mean_data, length(train_times), 1);
test_sx_all = test_sx_all - repmat(mean_data, length(test_times), 1);

%% PCA
[V, D] = PCA(train_sx_all, 2);

train_coeff = train_sx_all * V;
test_coeff = test_sx_all * V;

%% plot PCA projection
figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
idx = find(test_times > 0);
plot(test_coeff(idx,1),test_coeff(idx,2),'.')

figure;
scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
hold on
plot(test_coeff(:,1),test_coeff(:,2),'.')

return

%%
% figure;
% for i=1:length(test_times)
%     subplot(1,2,1)
%     cla
%     scatter(train_coeff(:,1), train_coeff(:, 2), 50, train_times)
%     hold on
%     plot(test_coeff(i,1),test_coeff(i,2),'.k','markersize', 20)
%
%     subplot(1,2,2)
%     imshow(test_images(:,:,i))
%     pause
%
% end

%% DMAPS

W = squareform(pdist(train_coeff)).^2;
eps = median(W(:))/10;
[V_dmaps, D_dmaps] = dmaps(W, eps, 10);

figure;
plot(V_dmaps(:,2),train_times,'.')

corr(V_dmaps(:,2),train_times, 'type', 'spearman')

%%
W = squareform(pdist([train_coeff; test_coeff])).^2;
eps = median(W(:))/30;
[V_dmaps, D_dmaps] = dmaps(W, eps, 10);

figure;
plot(V_dmaps(:,2),[train_times; test_times],'.')


%%

Wuu = squareform(pdist(test_coeff)).^2;
Wul = pdist2(test_coeff, train_coeff).^2;

eps = median(Wuu(:))/20;
Wuu = exp(-Wuu/eps);
Wul = exp(-Wul/eps);

fl = train_times - mean(train_times);
Du = diag(sum(Wuu));

fu = inv(Du - Wuu) * Wul * fl;

figure;
plot(test_times, fu, '.')




