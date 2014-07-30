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

for i=1:nmovies
    figure;
    subplot_dim1 = ceil(sqrt(nimages(i)));
    subplot_dim2 = ceil(nimages(i) / subplot_dim1);
    
    for j=1:nimages(i)
        subplot(subplot_dim1, subplot_dim2, j)
        imshow(images{i}(:,:,j));
    end
    
end

%%

PCA_data = cell(nmovies, 1);
PCA_data_mean = zeros(nmovies, npixels^2);

for i=1:nmovies
    PCA_data_tmp = reshape(double(images{i}), npixels^2, [])';
    PCA_data_mean(i, :) = mean(PCA_data_tmp);
    %     PCA_data_tmp = PCA_data_tmp - repmat(PCA_data_mean(i, :), nimages(i), 1);
    
    PCA_data{i} = PCA_data_tmp;
end

%%

err = zeros(nmovies);

for train_movie=1:nmovies
    
    train_movie  = 1
    nmodes = 10;
    
    h1 = figure;
    
    for i=1:nmovies
        pred_time = predict_times_PCA(PCA_data{train_movie}, time{train_movie}, PCA_data{i}, nmodes);
        
        err(train_movie, i) = mean(abs(pred_time - time{i}));
        
        figure(h1)
        subplot(3, 2, i)
        plot(time{i}, pred_time, '.')
        hold on
        plot(time{i}, time{i}, '-r')
        xlabel('true time')
        ylabel('predicted time')
        
    end
end

figure;
imagesc(err)
colorbar
xlabel('predicted movie')
ylabel('training movie')
title('average error')

%%

fixed_image_file = '../image_analysis_paper/late_images_aligned.mat';

load(fixed_image_file);

nfixed_images = size(image_set_raw_aligned_ordered, 4);
fixed_images = zeros(npixels, npixels, nfixed_images, 'uint8');
fixed_image_angle = 0;

for i=1:nfixed_images
    im_tmp = imresize(image_set_raw_aligned_ordered(:, :, channel, i), [npixels npixels]);
    im_tmp = image_fn(im_tmp, fixed_image_angle);
    
    fixed_images(:, :, i) = im_tmp;
end

figure;
subplot_dim1 = ceil(sqrt(nfixed_images));
subplot_dim2 = ceil(nfixed_images / subplot_dim1);
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
end

%%

train_movie = 1;

PCA_fixed_data = reshape(double(fixed_images), npixels^2, [])';

pred_time = predict_times_PCA(PCA_data{train_movie}, time{train_movie}, PCA_fixed_data, nmodes);

figure;
plot(pred_time, '.')
xlabel('image index (from diffusion maps)')
ylabel('predicted time')
title('Estimating times of fixed images')

[~, I] = sort(pred_time);
figure;
subplot_dim1 = ceil(sqrt(nfixed_images));
subplot_dim2 = ceil(nfixed_images / subplot_dim1);
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(image_set_raw_aligned_ordered(:,:,:,I(i)));
end


