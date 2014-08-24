clear all
close all

%% read in images

fixed_image_dir = '../../membrane_pictures/composite_108_noyolk';
nfixed_images = 108;

movies = {'../bomyi_emb01_gast01.avi';...
    '../bomyi_emb02_gast02.avi';...
    '../14_0623/emb01_hisRFP_gastrulation.avi';...
    '../14_0623/emb02_hisRFP_gastrulation.avi';...
    '../14_0624/emb01_hisRFP_gastrulation.avi';...
    '../14_0624/emb02_hisRFP_gastrulation.avi'; ...
    '../0709_emb02_cell_gast.avi'};

movies_angle = [-95;
    -55;
    -80;
    -80;
    -95;
    -120;
    -95];

movies_start = [20 25 15 10 15 5 110];
movies_end = [55 60 48 42 42 40 140];
dt = 0.5;
movie_idx = 3;

npixels = 500;
channel = 1;
nchannels = 3;

fixed_images = zeros(npixels, npixels, nchannels, nfixed_images, 'uint8');
fixed_image_angle = 0;

for i=1:nfixed_images
    im_tmp = imread(sprintf('%s/ordered%02d.tif', fixed_image_dir, i));
    
    fixed_images(:, :, :, i) = imresize(im_tmp, [npixels npixels]);
end
i = setdiff(1:nfixed_images, [79 92 94 96 99 100 102:108]);

fixed_images = fixed_images(:,:,:,i);
nfixed_images = length(i);

[movie_images, movie_times] = read_video(movies{movie_idx}, npixels);
movie_images = movie_images(:,:,:,movies_start(movie_idx):movies_end(movie_idx));
movie_times = dt * movie_times(movies_start(movie_idx):movies_end(movie_idx));
nmovie_images = length(movie_times);

movie_times = movie_times - 8;

for i=1:nmovie_images
    movie_images(:,:,:,i) = imrotate(movie_images(:,:,:,i), movies_angle(movie_idx), 'crop');
end

%%

movie_times_est = zeros(size(movie_times));
for i=1:nmovie_images
    movie_times_est(i) = predict_time_image(movie_images(:,:,channel, i));
end

figure;
plot(movie_times, movie_times_est,'.')
hold on 
plot([0 10], [0 10], '-r')

%% predict times

pred_times = zeros(nfixed_images, 1);
pred_times_min = zeros(nfixed_images, 1);
pred_times_max = zeros(nfixed_images, 1);

for i=1:nfixed_images
    [pred_times(i), pred_times_min(i), pred_times_max(i)] = predict_time_image(fixed_images(:,:,channel, i));
end

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
    
