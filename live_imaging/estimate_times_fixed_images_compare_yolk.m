clear all
close all

%%

channel = 1;
npixels = 100;

dt = 0.5;

crop_threshold = 1e3;

image_equalize = @(image, rot_angle) crop_image(adapthisteq(imrotate(image, rot_angle, 'crop')), crop_threshold, 0);
image_symmetrize = @(image) imlincomb(0.5, image, 0.5, fliplr(image));
image_fn = @(image, rot_angle) adapthisteq(image_symmetrize(image_equalize(image, rot_angle)));

movies = {'bomyi_emb01_gast01.avi';...
    'bomyi_emb02_gast02.avi';...
    '14_0623/emb01_hisRFP_gastrulation.avi';...
    '14_0623/emb02_hisRFP_gastrulation.avi';...
    '14_0624/emb01_hisRFP_gastrulation.avi';...
    '14_0624/emb02_hisRFP_gastrulation.avi'};

nmovies = length(movies);

theta = [-95;
    -55;
    -80;
    -80;
    -95;
    -120];

images = cell(nmovies, 1);
time = cell(nmovies, 1);
nimages = zeros(nmovies, 1);

%%
figure;
for i=1:nmovies
    [images_tmp, time{i}] = read_video(movies{i}, npixels);
    
    time{i} = time{i} * dt;
    
    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);
    
    nimages(i) = length(time{i});
    
    for j=1:nimages(i)
        images_tmp(:, :, j) = image_fn(images_tmp(:, :, j), theta(i));
    end
    
    images{i} = images_tmp;
    
    subplot(2, 3, i)
    imshow(images{i}(:,:, round(nimages(i)/2)));
end

%%

PCA_data = cell(nmovies, 1);
for i=1:nmovies
    PCA_data{i} = create_PCA_data(images{i});
end

%%

nmodes = 8;

shift_times = estimate_shift_times(PCA_data, time, nmodes);

time_adjust = time;
for i=1:nmovies
    time_adjust{i} = time_adjust{i} - shift_times(i);
end

%%

fixed_image_dir = '../membrane_pictures/composite_108_noyolk';
nfixed_images = 108;
subplot_dim1 = 6;
subplot_dim2 = 18;

fixed_images = zeros(npixels, npixels, nfixed_images, 'uint8');
fixed_image_angle = 0;

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    im_tmp = imresize(im_tmp(:, :, channel), [npixels npixels]);
    im_tmp = image_fn(im_tmp, fixed_image_angle);
    
    fixed_images(:, :, i) = im_tmp;
end

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
end

%%

PCA_fixed_data = create_PCA_data(fixed_images);

PCA_data_movies = vertcat(PCA_data{:});
time_movies = vertcat(time_adjust{:});

pred_time_all = predict_times_PCA(PCA_data_movies, time_movies, PCA_fixed_data, nmodes);

figure;
plot(pred_time_all, '.')
xlabel('image index (from diffusion maps)')
ylabel('predicted time')
title('Estimating times of fixed images')

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
    axis tight
    title(sprintf('t = %2.2f', pred_time_all(i)))
end

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    imshow(im_tmp);
    axis tight
    title(sprintf('t = %2.2f', pred_time_all(i)))
end

%%

time_monotonic = fit_monotonic_curve(pred_time_all);

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
    axis tight
    title(sprintf('t = %2.2f', time_monotonic(i)))
end

%%

fixed_image_file2 = '../image_analysis_paper/late_images_aligned.mat';

fixed_images2 = zeros(npixels, npixels, nfixed_images, 'uint8');


load(fixed_image_file2);

for i=1:nfixed_images
    im_tmp = image_set_raw_aligned_ordered(:, :, :, i);
    im_tmp = imresize(im_tmp(:, :, channel), [npixels npixels]);
    im_tmp = image_fn(im_tmp, fixed_image_angle);
    
    fixed_images2(:, :, i) = im_tmp;
end

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images2(:,:,i));
end

%%

PCA_fixed_data2 = create_PCA_data(fixed_images2);

pred_time_all2 = predict_times_PCA(PCA_data_movies, time_movies, PCA_fixed_data2, nmodes);

figure;
plot(pred_time_all2, '.')
xlabel('image index (from diffusion maps)')
ylabel('predicted time')
title('Estimating times of fixed images')

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images2(:,:,i));
    axis tight
    title(sprintf('t = %2.2f', pred_time_all2(i)))
end

figure; 
plot(pred_time_all, pred_time_all2, '.')