clear all
close all

%% read in movie
file_name = {'bomyi_emb01_gast01.avi'; 'bomyi_emb02_gast02.avi'};
num_files = 2;

npixels = 100;
nchannels = 1;
channel = 1;

times = [];
num_images = 0;
num_images_vec = [];
image_set = [];

for j=1:num_files
    
    vidobj = VideoReader(file_name{j});

    num_images_tmp = vidobj.NumberOfFrames;

    image_set_tmp = zeros(npixels, npixels, num_images_tmp, 'uint8');
%     image_set_tmp = true(npixels, npixels, num_images_tmp);
    
    for i = 1:num_images_tmp
        A = read(vidobj, i);
        A =  A(:, :, channel);
%             A = im2bw(A, 0.15);
        A = imresize(A, [npixels npixels]);
        image_set_tmp(:, :, i) = A;
    end
    
    image_set = cat(3, image_set, image_set_tmp);
    times = [times; (1:num_images_tmp)'];
    num_images = num_images + num_images_tmp;
    num_images_vec = [num_images_vec; num_images_tmp];
end

subplot_dim1 = floor(sqrt(num_images));
subplot_dim2 = ceil(num_images/subplot_dim1);


%%

% figure;
% for i=1:num_images
%     imshow(image_set(:, :, i));
%     pause(0.01);
%     clf
% end


%%
figure;
for i=1:num_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(image_set(:, :, i));
end


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
x = double(image_set(:,:,1));
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
sx_all = zeros(num_images, size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:num_images
    i
    x = double(image_set(:,:,i));
    
    Sx = scat(x, Wop);
    
    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%%

W = squareform(pdist(sx_all)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

figure;
plot(V(:,2),times, '.');


%%
[~, I] = sort(V(:,2));
figure;
for i=1:num_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(image_set(:, :, I(i)));
end

%%

figure;
scatter(V(:,2),V(:,4),50, times, '.')
