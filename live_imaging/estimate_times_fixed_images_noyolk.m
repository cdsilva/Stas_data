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

% for i=1:nmovies
%     figure;
%     subplot_dim1 = ceil(sqrt(nimages(i)));
%     subplot_dim2 = ceil(nimages(i) / subplot_dim1);
%
%     for j=1:nimages(i)
%         subplot(subplot_dim1, subplot_dim2, j)
%         imshow(images{i}(:,:,j));
%     end
%
% end

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

bias3 = zeros(nmovies, 1);
variance3 = zeros(nmovies, 1);

% figure;
for i=1:nmovies
    
    data_tmp = vertcat(PCA_data{1:nmovies ~= i});
    time_tmp = vertcat(time_adjust{1:nmovies ~= i});
    
    pred_time = predict_times_PCA(data_tmp, time_tmp, PCA_data{i}, nmodes);
    
    bias3(i) = mean(pred_time - time_adjust{i});
    variance3(i) = mean((pred_time - time_adjust{i}).^2);
    
    %     subplot(3, 2, i)
    %     plot(time_adjust{i}, pred_time, '.')
    %     hold on
    %     plot(time_adjust{i}, time_adjust{i}, '-r')
    %     xlabel('true time')
    %     ylabel('predicted time')
    
end

figure;
bar(bias3)

figure;
bar(variance3.^0.5)

fprintf('cross validation average error : %2.2f \n', sqrt(mean(variance3)));

%%

fixed_image_dir = '../membrane_pictures/composite_50_noyolk';
nfixed_images = 50;
fixed_images = zeros(npixels, npixels, nfixed_images, 'uint8');
fixed_image_angle = 0;


for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/emb%02d.tif', fixed_image_dir, i));
    im_tmp = imresize(im_tmp(:, :, channel), [npixels npixels]);
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
% figure;
%
% pred_time_fixed_images = zeros(nfixed_images, nmovies);
% PCA_fixed_data = reshape(double(fixed_images), npixels^2, [])';
%
% for i = 1:nmovies;
%
%
%     pred_time(:, i) = predict_times_PCA(PCA_data{i}, time_adjust{i}, PCA_fixed_data, nmodes);
%
%     subplot(3,2,i)
%     plot(pred_time(:, i), '.')
%     xlabel('image index (from diffusion maps)')
%     ylabel('predicted time')
%     title('Estimating times of fixed images')
%
% end
%
% figure;
% for i=1:nmovies
%     for j=1:nmovies
%         subplot(nmovies, nmovies, (i-1)*nmovies+j)
%         plot(pred_time(:, i), pred_time(:,j), '.')
%         hold on
%         plot(pred_time(:, i), pred_time(:,i), '-r')
%     end
% end

%%

PCA_fixed_data = create_PCA_data(fixed_images);

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

figure;
subplot_dim1 = 5;
subplot_dim2 = 10;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(fixed_images(:,:,i));
    axis tight
    title(sprintf('t = %2.2f', pred_time_all(i)))
end

figure;
for i=1:nfixed_images
    subplot(subplot_dim1, subplot_dim2, i)
    im_tmp = imread(sprintf('%s/emb%02d.tif', fixed_image_dir, i));    
    imshow(im_tmp);
    axis tight
    title(sprintf('t = %2.2f', pred_time_all(i)))
end

