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
scatnet_path = 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
% scatnet_path = '../../../MATLAB/scatnet-0.2';

% directory name where images are stored
dir_name = '../membrane_pictures/14_0403_dpERK_late';

nimages = 46;

%% define function to adjust each image
image_fn = @(im1) adapthisteq(medfilt2(imresize(im1, [100 100]), [5 5]),'distribution','exponential');


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
% sx_all = zeros(nimages, size(Sx_mat, 1), nchannels);

% % compute scattering invariants for each image
% for i=1:nimages
%     i
%     for j=1:nchannels
%         x = double(fixed_images(:,:,j,i));
%
%         Sx = scat(x, Wop);
%
%         sx_all(i, :, j) = mean(mean(format_scat(Sx),2),3)';
%     end
% end

sx_all = zeros(nimages, size(Sx_mat, 1));
% compute scattering invariants for each image
for i=1:nimages
    i
    x = double(fixed_images(:,:,1,i));
    for j=2:nchannels
        x = x + 2*double(fixed_images(:,:,j,i));
    end
    
    Sx = scat(x, Wop);
    
    sx_all(i, :) = mean(mean(format_scat(Sx),2),3)';
end

%%

W = squareform(pdist(sx_all)).^2;

% W = zeros(nimages);
% weights = [1 0.5 0.5];
%
% for i=1:nchannels
%     W = W + weights(i)*squareform(pdist(sx_all(:,:,i))).^2;
% end

%%

eps = median(W(:))/10;
[V, D] = dmaps(W, eps, 10);

[~, I] = sort(V(:,2));

figure;
for i=1:nimages
    subplot(7,7,i)
    imshow(fixed_images(:,:,:,I(i)))
end
