function [images, time] = load_fixed_images1(npixels)
% images with nuclei (red), dpERK (green), and Dl (blue) during
% gastrulation

fixed_image_dir = '../membrane_pictures/composite_108_noyolk';
nfixed_images = 108;

nchannels = 3;
channel = 1;

images = zeros(npixels, npixels, nchannels, nfixed_images, 'uint8');
time = zeros(nfixed_images, 1);

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    
    [time(i), ~, ~] = predict_time_image(imresize(im_tmp(:,:,channel), [100 100]));
    
    images(:,:,:,i) = crop_hires_image(imresize(im_tmp, [npixels npixels]), channel);
end
idx = setdiff(1:nfixed_images, [79 92 94 96 99 100 102:108]);

images = images(:,:,:,idx);
time = time(idx);