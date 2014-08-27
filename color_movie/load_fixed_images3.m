function [images, time] = load_fixed_images3(npixels)
% images with dpERK (red), ind (green), and nuclei (blue) during
% cellularization and gastrulation

nchannels = 3;
channel = 3;

fixed_image_dir = '../membrane_pictures/ind_carmeline';
nfixed_images = 140;

load(sprintf('%s/time.mat', fixed_image_dir));
time = time_all;

images = zeros(npixels, npixels, nchannels, nfixed_images, 'uint8');

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/emb%02d.tif', fixed_image_dir, i));
    images(:,:,:,i) = crop_hires_image(imresize(im_tmp, [npixels npixels]), channel);
end

idx = setdiff(1:nfixed_images, [17 49 76 98 121 124 129 131 136 139 140]);
images = images(:,:,:,idx);
time = time(idx);
