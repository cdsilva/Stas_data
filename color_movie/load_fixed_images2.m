function [images, time] = load_fixed_images2(npixels)
% images with nuclei (red), dpERK (green), and Dl (blue) during
% cellularization

nchannels = 3;
channel = 1;

fixed_image_dir = '../membrane_pictures/large_dataset';
nfixed_images = 90;

% load('calibration_curve_data.mat');
load(sprintf('%s/mem_lengths.mat', fixed_image_dir));
% clear length;
% time = interp1(cal_curve_membrane, cal_curve_time, mem_lengths);
time = mem_lengths;

npixels2 = 100;

images = zeros(npixels, npixels, nchannels, nfixed_images, 'uint8');
images_tmp = zeros(npixels2, npixels2, nchannels, nfixed_images, 'uint8');

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/lat%02d.tif', fixed_image_dir, i));
    images(:,:,:,i) = imresize(im_tmp, [npixels npixels]);
    images_tmp(:,:,:,i) = crop_image(imresize(im_tmp, [npixels2 npixels2]), channel);
end

images_tmp(:,:,1,:) = 0;

nrot = 10;
nshifts = 0;
shift_max = 0;
display_waitbar = true;

[R, W] = compute_pairwise_alignments(images_tmp, nrot, nshifts, shift_max, display_waitbar);
dim = size(R, 1) / nfixed_images;

%
neigs = 10;
alpha = 0;
W2 = W.^2;
eps = median(W2(:))/10;
[R_opt, ~, ~, ~] = vdm(R, W2, eps, neigs, alpha);

%
for i=1:nfixed_images
    R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
    theta_opt = atan2d(R_tmp(2,1), R_tmp(1,1));
    images(:,:,:,i) = crop_hires_image(imrotate(images(:,:,:,i), -theta_opt-45, 'crop'), channel);
end

idx = setdiff(1:nfixed_images, [5 7 8 10 18 20 22 36 65 70]);
images = images(:,:,:,idx);
time = time(idx);


