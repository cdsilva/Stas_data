clear all
close all

%% read in images

fixed_image_dir = '../../membrane_pictures/composite_108_noyolk';
nfixed_images = 108;

fixed_image_dir2 = '../../membrane_pictures/large_dataset';
nfixed_images2 = 90;

% movie_file = '../bomyi_emb01_gast01.avi';
% movie_angle = -95;
% movie_start = 20;
% movie_end = 55;
% 
% movie_file = '../bomyi_emb02_gast02.avi';
% movie_angle = -55;
% movie_start = 25;
% movie_end = 60;
% 
movie_file = '../14_0623/emb01_hisRFP_gastrulation.avi';
movie_angle = -80;
movie_start = 1;
movie_end = 48;
% 
% movie_file = '../14_0623/emb02_hisRFP_gastrulation.avi';
% movie_angle = -80;
% movie_start = 10;
% movie_end = 42;
% 
% movie_file = '../14_0624/emb01_hisRFP_gastrulation.avi';
% movie_angle = -95;
% movie_start = 15;
% movie_end = 42;
% 
% movie_file = '../14_0624/emb02_hisRFP_gastrulation.avi';
% movie_angle = -120;
% movie_start = 5;
% movie_end = 40;
% 
% movie_file = '../0709_emb02_cell_gast.avi';
% movie_angle = -95;
% movie_start = 110;
% movie_end = 140;

dt = 0.5;

npixels = 500;
channel = 1;
nchannels = 3;


%%
[movie_images, movie_times] = read_video(movie_file, npixels);
movie_images = movie_images(:,:,:,movie_start:movie_end);
movie_times = dt * movie_times(movie_start:movie_end);
nmovie_images = length(movie_times);

for i=1:nmovie_images
    movie_images(:,:,:,i) = imrotate(movie_images(:,:,:,i), movie_angle, 'crop');
end

%%

movie_times_est = zeros(size(movie_times));
for i=1:nmovie_images
    movie_times_est(i) = predict_time_image(movie_images(:,:,channel, i));
end

time_est_offset = 10;
movie_times_est = movie_times_est + time_est_offset;

figure;
plot(movie_times, movie_times_est,'.')
hold on
plot([0 10], [0 10], '-r')

return

%%
fixed_images = [];
pred_times = [];

%%
fixed_image_angle = 0;

fixed_images_tmp = zeros(npixels, npixels, nchannels, nfixed_images, 'uint8');
for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    
    fixed_images_tmp(:, :, :, i) = imresize(im_tmp, [npixels npixels]);
end
i = setdiff(1:nfixed_images, [79 92 94 96 99 100 102:108]);

fixed_images_tmp = fixed_images_tmp(:,:,:,i);
nfixed_images = length(i);

fixed_images = cat(4, fixed_images, fixed_images_tmp);

% predict times

pred_times_tmp = zeros(nfixed_images, 1);

for i=1:nfixed_images
    [pred_times_tmp(i), ~, ~] = predict_time_image(fixed_images_tmp(:,:,channel, i));
end

pred_times = [pred_times; pred_times_tmp];

%%

fixed_images_tmp = zeros(npixels, npixels, nchannels, nfixed_images2, 'uint8');
for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/lat%02d.tif', fixed_image_dir2, i));
    
    fixed_images_tmp(:, :, :, i) = imresize(im_tmp, [npixels npixels]);
end

fixed_images = cat(4, fixed_images, fixed_images_tmp);

% predict times

pred_times_tmp = zeros(nfixed_images, 1);

for i=1:nfixed_images
    [pred_times_tmp(i), ~, ~] = predict_time_image(fixed_images_tmp(:,:,channel, i));
end

pred_times = [pred_times; pred_times_tmp];

%% plot
figure;
plot(pred_times,'.')


%%
fixed_images_cropped = fixed_images;
for i=1:nfixed_images
    fixed_images_cropped(:,:,:,i) = crop_hires_image(fixed_images(:,:,:,i), channel);
end

movie_images_colored = movie_images;
kernel_scale = 0.5;

for i=1:nmovie_images
    weights = exp(-(pred_times - movie_times(i)).^2/kernel_scale.^2);
    weights = weights / sum(weights);
    
    movie_images_colored(:,:,:,i) = immultiply(fixed_images_cropped(:,:,:,1), weights(1));
    for j=2:nfixed_images
        movie_images_colored(:,:,:,i) = imlincomb(1, movie_images_colored(:,:,:,i), weights(j), fixed_images_cropped(:,:,:,j));
    end
    movie_images_colored(:,:,channel,i) = adapthisteq(crop_hires_image(movie_images(:,:,channel, i), 0), 'numtiles', [16 16]);
end

%%

writerObj = VideoWriter('colored_movie.avi');
writerObj.FrameRate = 10;
open(writerObj);

figure;
% set(gcf, 'paperunits', 'centimeters')
% set(gcf, 'papersize', [8 8])
% set(gcf, 'paperposition',[0 0 8 8])

imshow(make_gray_nuclei(movie_images_colored(:,:,:,1)), 'initialmagnification','fit','border','tight')
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
clf

for i=1:nmovie_images
    imshow(make_gray_nuclei(movie_images_colored(:,:,:,i)))
    
    frame = getframe;
    writeVideo(writerObj,frame);
    clf
end

close(writerObj);

