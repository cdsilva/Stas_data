clear all
close all

%% read in images

fixed_image_dir = '../../membrane_pictures/composite_108_noyolk';
nfixed_images = 108;
npixels = 100;
channel = 1;
fixed_images = zeros(npixels, npixels, nfixed_images, 'uint8');
fixed_image_angle = 0;

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    
    fixed_images(:, :, i) = imresize(im_tmp(:,:,channel), [npixels npixels]);
end

figure;
subplot_dim1 = ceil(sqrt(nfixed_images));
subplot_dim2 = ceil(nfixed_images / subplot_dim1);
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
end

%% predict times

pred_times = zeros(nfixed_images, 1);
pred_times_min = zeros(nfixed_images, 1);
pred_times_max = zeros(nfixed_images, 1);

for i=1:nfixed_images
    [pred_times(i), pred_times_min(i), pred_times_max(i)] = predict_time_image(fixed_images(:,:,i));
end

%% plot
figure;
plot(pred_times,'.')

figure;
plot(pred_times_min)
hold on
plot(pred_times_max, '-r')

figure;
plot(pred_times_max-pred_times_min)



