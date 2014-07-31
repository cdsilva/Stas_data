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

bias = zeros(nmovies);
variance = zeros(nmovies);

for train_movie=1:nmovies
    
    nmodes = 10;
    
    h1 = figure;
    
    for i=1:nmovies
        pred_time = predict_times_PCA(PCA_data{train_movie}, time{train_movie}, PCA_data{i}, nmodes);
        
        bias(train_movie, i) = mean(pred_time - time{i});
        variance(train_movie, i) = mean((pred_time - time{i}).^2);
        
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
imagesc(bias)
colorbar
xlabel('predicted movie')
ylabel('training movie')
title('average bias')

figure;
imagesc(variance.^0.5)
colorbar
xlabel('predicted movie')
ylabel('training movie')
title('average std dev')

%% adjust movies

synch_matrix = (bias-bias')/2;
synch_scale = 4*(max(synch_matrix(:))-min(synch_matrix(:)));
synch_matrix = exp(sqrt(-1) * synch_matrix * (2*pi/synch_scale));
[V, D] = eigs(synch_matrix, 1);
shift_times = atan2(imag(V), real(V)) * synch_scale/(2*pi);
shift_times = shift_times - max(shift_times);

time_adjust = time;
for i=1:nmovies
    %     time_adjust{i} = time_adjust{i} + mean(bias(1:nmovies ~= i, i));
    time_adjust{i} = time_adjust{i} - shift_times(i);
end

%% retrain

bias2 = zeros(nmovies);
variance2 = zeros(nmovies);

for train_movie=1:nmovies
    
    nmodes = 10;
    
    h1 = figure;
    
    for i=1:nmovies
        pred_time = predict_times_PCA(PCA_data{train_movie}, time_adjust{train_movie}, PCA_data{i}, nmodes);
        
        bias2(train_movie, i) = mean(pred_time - time_adjust{i});
        variance2(train_movie, i) = mean((pred_time - time_adjust{i}).^2);
        
        figure(h1)
        subplot(3, 2, i)
        plot(time_adjust{i}, pred_time, '.')
        hold on
        plot(time_adjust{i}, time_adjust{i}, '-r')
        xlabel('true time')
        ylabel('predicted time')
        
    end
end

figure;
imagesc(bias2)
colorbar
xlabel('predicted movie')
ylabel('training movie')
title('average bias')

figure;
imagesc(variance2.^0.5)
colorbar
xlabel('predicted movie')
ylabel('training movie')
title('average std dev')

%%

bias3 = zeros(nmovies, 1);
variance3 = zeros(nmovies, 1);

figure;
for i=1:nmovies
    
    nmodes = 10;
    
    data_tmp = [];
    time_tmp = [];
    for j=1:nmovies
        if i ~= j
            data_tmp = [data_tmp; PCA_data{j}];
            time_tmp = [time_tmp; time_adjust{j}];
        end
    end
    
    pred_time = predict_times_PCA(data_tmp, time_tmp, PCA_data{i}, nmodes);
    
    bias3(i) = mean(pred_time - time_adjust{i});
    variance3(i) = mean((pred_time - time_adjust{i}).^2);
    
    subplot(3, 2, i)
    plot(time_adjust{i}, pred_time, '.')
    hold on
    plot(time_adjust{i}, time_adjust{i}, '-r')
    xlabel('true time')
    ylabel('predicted time')
    
end

figure;
bar(bias3)

figure;
bar(variance3.^0.5)

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
figure;

pred_time_fixed_images = zeros(nfixed_images, nmovies);
PCA_fixed_data = reshape(double(fixed_images), npixels^2, [])';

for i = 1:nmovies;
    
    
    pred_time(:, i) = predict_times_PCA(PCA_data{i}, time_adjust{i}, PCA_fixed_data, nmodes);
    
    subplot(3,2,i)
    plot(pred_time(:, i), '.')
    xlabel('image index (from diffusion maps)')
    ylabel('predicted time')
    title('Estimating times of fixed images')
    
end

figure;
for i=1:nmovies
    for j=1:nmovies
        subplot(nmovies, nmovies, (i-1)*nmovies+j)
        plot(pred_time(:, i), pred_time(:,j), '.')
        hold on
        plot(pred_time(:, i), pred_time(:,i), '-r')
    end
end

%%

data_tmp = [];
time_tmp = [];
for j=1:nmovies
    data_tmp = [data_tmp; PCA_data{j}];
    time_tmp = [time_tmp; time_adjust{j}];
end

pred_time_all = predict_times_PCA(data_tmp, time_tmp, PCA_fixed_data, nmodes);

figure;
plot(pred_time_all, '.')
xlabel('image index (from diffusion maps)')
ylabel('predicted time')
title('Estimating times of fixed images')
