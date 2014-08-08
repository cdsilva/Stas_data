clear all
% close all

%%

channel = 1;

fixed_image_dir = '../membrane_pictures/composite_108_noyolk';
nfixed_images = 108;
subplot_dim1 = 6;
subplot_dim2 = 18;

bomyi_idx = [6 4 5 2 1 18 3 8 9 16 13 12 7 11 14 15 17 19 20 22 ...
    29 38 32 27 24 25 10 21 23 26 28 30 31 33 34 36 35 37 41 40 ...
    39 43 42 45 49 53 55 46 48 54 52 44 47 50 56 61 51 59 65 78 ...
    72 66 58 60 70 81 64 57 67 95 98 77 76 83 90 87 82 74 73 62 ...
    63 80 91 68 71 88 84 93 92 94 101 85 97 89 86 75 100 107 69 103 ...
    108 79 104 106 99 96 102 105];

im_tmp = imread(sprintf('%s/ordered01.tif', fixed_image_dir));
npixels = size(im_tmp, 1);
fixed_images = zeros(npixels, npixels, nfixed_images, 'uint8');

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    
    fixed_images(:, :, i) = im_tmp(:, :, channel);
end

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
end

%%

t_pred = zeros(nfixed_images, 1);
tmin_pred = zeros(nfixed_images, 1);
tmax_pred = zeros(nfixed_images, 1);

for i=1:nfixed_images
    [t_pred(i), tmin_pred(i), tmax_pred(i)] = predict_time_image(fixed_images(:,:,i));
end

%%

figure;
errorbar(1:nfixed_images, t_pred(bomyi_idx), t_pred(bomyi_idx)-tmin_pred(bomyi_idx), tmax_pred(bomyi_idx)-t_pred(bomyi_idx), '.b')

