clear all
close all

%% read in movies

npixels = 100;
channel = 1;

file_name = '14_0624/emb01_cell_gast.avi';

[images, time] = read_video(file_name, npixels);
images = images(:, :, channel, :);
images = squeeze(images);

%% adjust images

image_fn = @(image) adapthisteq(medfilt2(imresize(image, [100 100]), [5 5]),'distribution','exponential');

%% compute scattering coefficients

addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
% addpath '../../../MATLAB/scatnet-0.2';
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

%% mean-center
mean_data = mean(sx_all);

sx_all = sx_all - repmat(mean_data, length(time), 1);

%% PCA
[V, D] = PCA(sx_all, 100);

if mean(V(:,1)) > 0
    V(:,1) = -V(:,1);
end

figure;
bar(diag(D(1:10,1:10))/sum(diag(D)))

coeff = sx_all * V;

%% plot PCA projection

figure;
scatter(coeff(:,1), coeff(:, 2), 50, time, '.')

figure;
scatter3(coeff(:,1), coeff(:, 2),coeff(:, 3), 50, time, '.')





