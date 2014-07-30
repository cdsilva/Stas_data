clear all
close all


%% parameters

% number of pixels to subsample the movies to (doing the scattering
% transform on the original high-resolution movies takes too long)
npixels = 100;

% channel of the movie which contains the revevant signal (red = 1, green =
% 2, blue = 3)
channel = 1;

% file names of movies
file_name1 = 'bomyi_emb01_gast01.avi';
file_name2 = 'bomyi_emb02_gast02.avi';

image_fn = @(image) medfilt2(adapthisteq(crop_image(image, 100, 0)),[5 5]);

%% read in movies

[images1, times1] = read_video(file_name1, npixels);
images1 = images1(:, :, channel, :);
images1 = squeeze(images1);

[images2, times2] = read_video(file_name2, npixels);
images2 = images2(:, :, channel, :);
images2 = squeeze(images2);

for i=1:length(times1)
    images1(:,:,i) = image_fn(images1(:,:,i));
    images1(:,:,i) = adapthisteq(images1(:,:,i)-uint8(0.75*mean(double(images1), 3)));
end

for i=1:length(times2)
    images2(:,:,i) = image_fn(images2(:,:,i));
    images2(:,:,i) = adapthisteq(images2(:,:,i)-uint8(0.75*mean(double(images1), 3)));
end

%%

figure;
for i=1:length(times1)
    subplot(10,10,i)
    imshow(images1(:,:,i)); 
end

figure;
for i=1:length(times2)
    subplot(10,10,i)
    imshow(images2(:,:,i));
end

%% compare frames of movies

test_idx = 50;

figure;
subplot(1,2,1)
imshow(images1(:,:,test_idx))
subplot(1,2,2)
imshow(images2(:,:,test_idx))

%% rotate and shift movies to consistent frame of reference

theta = 45;
shift = [-15 -9];

for i=1:length(times2)
    images2(:,:,i) = imrotate(images2(:,:,i), theta, 'crop');
%     images2(:,:,i) = circshift(images2(:,:,i), shift);
end

%% compare frames after rotation and shift

im_tmp = zeros(npixels, npixels, 3, 'uint8');
im_tmp(:,:,1) = images1(:,:,test_idx);
im_tmp(:,:,3) = images2(:,:,test_idx);

figure;
imshow(im_tmp)
title('images from two movies overlaid to compare alignment')

%% compute pairwise distances

data1 = reshape(double(images1), npixels^2, length(times1))';
data2 = reshape(double(images2), npixels^2, length(times2))';
pdist_matrix = pdist2(data1(randperm(length(times1)), :), data2);
% pdist_matrix = pdist2(data1, data2);

figure;
imagesc(pdist_matrix)
title('distance matrix between movies')

%% time warp

[D, warp_path] = DTW(pdist_matrix);

figure;
imagesc(D)
title('time warp matrix')

figure;
imagesc(warp_path)
colormap(gray)
title('time warp path')

D(end, end)

figure;
plot(sort(D(warp_path(:)>0)),'.')
