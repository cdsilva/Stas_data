clear all
close all

%%

channel = 1;
npixels = 100;

dt = 0.5;

crop_threshold = 1e3;

image_equalize = @(image, rot_angle) crop_image(adapthisteq(imrotate(image, rot_angle, 'crop'), 'cliplimit',0.005), crop_threshold, 0);
image_symmetrize = @(image) imlincomb(0.5, image, 0.5, fliplr(image));
image_fn = @(image, rot_angle) adapthisteq(image_symmetrize(image_equalize(image, rot_angle)));

% movies = {'bomyi_emb01_gast01.avi';...
%     'bomyi_emb02_gast02.avi';...
%     '14_0623/emb01_hisRFP_gastrulation.avi';...
%     '14_0623/emb02_hisRFP_gastrulation.avi';...
%     '14_0624/emb01_hisRFP_gastrulation.avi';...
%     '14_0624/emb02_hisRFP_gastrulation.avi'; ...
%     '0709_emb02_cell_gast.avi'};
%
% nmovies = length(movies);
%
% theta = [-95;
%     -55;
%     -80;
%     -80;
%     -95;
%     -120;
%     -95];
%
% movie_start = [1; 1; 1; 1; 1; 1; 96];
% movie_end = [20; 18; 8; 7; 0; 0; 0];
%
% make_subplot = @(i) subplot(2, 4, i);
%
% images = cell(nmovies, 1);
% time = cell(nmovies, 1);
% nimages = zeros(nmovies, 1);
%
% %%
% figure;
% for i=1:nmovies
%     [images_tmp, time_tmp] = read_video(movies{i}, npixels);
%     images_tmp = images_tmp(:, :, :, movie_start(i):end-movie_end(i));
%     time_tmp = time_tmp(movie_start(i):end-movie_end(i)) - time_tmp(movie_start(i));
%
%     time{i} = time_tmp * dt;
%
%     images_tmp = images_tmp(:, :, channel, :);
%     images_tmp = squeeze(images_tmp);
%
%     nimages(i) = length(time{i});
%
%     for j=1:nimages(i)
%         images_tmp(:, :, j) = image_fn(images_tmp(:, :, j), theta(i));
%     end
%
%     images{i} = images_tmp;
%
%     make_subplot(i);
%     imshow(images{i}(:,:, round(nimages(i)/2)));
% end
%
% % for i=1:nmovies
% %     figure;
% %     subplot_dim1 = ceil(sqrt(nimages(i)));
% %     subplot_dim2 = ceil(nimages(i) / subplot_dim1);
% %
% %     for j=1:nimages(i)
% %         subplot(subplot_dim1, subplot_dim2, j)
% %         imshow(images{i}(:,:,j));
% %     end
% %
% % end
%
% %%
%
% PCA_data = cell(nmovies, 1);
% for i=1:nmovies
%     PCA_data{i} = create_PCA_data(images{i});
% end
%
% %%
%
% nmodes = 5;
%
% shift_times = estimate_shift_times(PCA_data, time, nmodes);
%
% time_adjust = time;
% for i=1:nmovies
%     time_adjust{i} = time_adjust{i} - shift_times(i);
% end

%%

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
    imshow(imresize(fixed_images(:,:,i), [20 20]));
end

return
%%

npixels2 = 15;

[X, Y] = meshgrid(1:npixels2, 1:npixels2);
X = X(:);
Y = Y(:);
D = squareform(pdist([X Y]));

emd_dist = zeros(nfixed_images);
for i=1:nfixed_images
    i
    
    imagei = reshape(double(imresize(fixed_images(:, :, i), [npixels2 npixels2])), [], 1);
    for j=1:i-1
        
        imagej = reshape(double(imresize(fixed_images(:, :, j), [npixels2 npixels2])), [], 1);
        
        emd_dist(i, j) = emd_hat_gd_metric_mex(imagei, imagej, D);
    end
end

emd_dist = emd_dist + emd_dist';

figure;
imagesc(emd_dist(bomyi_idx, bomyi_idx))


