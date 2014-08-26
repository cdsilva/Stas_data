clear all
close all

%%

npixels = 100;
dt = 0.5;

% [images, time] = read_video('../live_imaging/bomyi_emb01_gast01.avi', npixels);
% save('movie1.mat', 'images','time');
% load('movie1.mat');
% idx = 21:56;

% [images, time] = read_video('../live_imaging/bomyi_emb02_gast02.avi', npixels);
% save('movie2.mat', 'images','time');
% load('movie2.mat');
% idx = 31:60;

% [images, time] = read_video('../live_imaging/14_0623/emb01_hisRFP_gastrulation.avi', npixels);
% save('movie3.mat', 'images','time');
% load('movie3.mat');
% idx = 27:47;

% [images, time] = read_video('../live_imaging/14_0623/emb02_hisRFP_gastrulation.avi', npixels);
% save('movie4.mat', 'images','time');
% load('movie4.mat');
% idx = 24:44;    
% FEW OUTLIERS

% [images, time] = read_video('../live_imaging/14_0624/emb01_hisRFP_gastrulation.avi', npixels);
% save('movie5.mat', 'images','time');
% load('movie5.mat');
% idx = 18:46;    

% [images, time] = read_video('../live_imaging/14_0624/emb02_hisRFP_gastrulation.avi', npixels);
% save('movie6.mat', 'images','time');
% load('movie6.mat');
% idx = 17:42;

% [images, time] = read_video('../live_imaging/0709_emb02_cell_gast.avi', 100);
% save('movie7.mat', 'images','time');

nimages = length(time);
make_subplot = @(i) subplot(ceil(nimages/10), 10, i);

image_set = zeros(npixels, npixels, nimages, 'uint8');
channel = 1;

time = time * dt;

filt = fspecial('disk');

for i=1:nimages
    image = images(:, :, channel, i);
    image = imresize(image, [npixels npixels]);

    image = crop_image(image, channel);
    image = padarray(image, [10 10]);
    image = imresize(image, [npixels npixels]);

    image = adapthisteq(image);

    image = imfilter(image, filt, 'replicate');
    
    image_set(:, :, i) = image;
end

%%

figure;
for i=1:nimages
    make_subplot(i)
    imshow(image_set(:, :, i))
end

%%

image_set = image_set(:, :, idx);
time = time(idx);
nimages = length(time);
make_subplot = @(i) subplot(ceil(nimages/10), 10, i);

%%
figure;
for i=1:nimages
    make_subplot(i)
    imshow(image_set(:, :, i))
end


%%

rng(123);
theta = 360*rand(nimages, 1);
% theta = randi(20, nimages, 1) * 360 / 20;
idx = randperm(nimages);

for i=1:nimages
    image_set(:, :, i) = imrotate(image_set(:, :, i), theta(i), 'crop');
end

image_set = image_set(:, :, idx);
time = time(idx);
theta = theta(idx);

figure;
for i=1:nimages
    make_subplot(i)
    imshow(image_set(:, :, i))
end

%%

nrot = 40;
nshifts = 0;
shift_max = 0.1;

tic
[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
toc

dim = size(R, 1) / nimages;

%%
W2 = W.^2;
eps = median(W2(:))/10;
neigs = 10;
alpha = 0;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);

%%
theta_opt = zeros(size(theta));
figure;
for i=1:nimages
    R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
    theta_opt(i) = atan2d(R_tmp(2,1), R_tmp(1,1));

    make_subplot(i)
    imshow(rotate_image(image_set(:, :, i), R_opt(dim*(i-1)+1:dim*i, :)'))
end

%%
figure;
plot(mod(theta-theta(1), 360), mod(theta_opt+10, 360)-10, '.')
xlabel('true angle')
ylabel('recovered angle')

%%
if corr(time, embed_coord(:,1)) < 0
    embed_coord(:,1) = -embed_coord(:,1);
end

[~, I] = sort(embed_coord(:, 1));

figure;
for i=1:nimages
    make_subplot(i)
    imshow(rotate_image(image_set(:, :, I(i)), R_opt(dim*(I(i)-1)+1:dim*I(i), :)'))
end

%%

figure;
plot(time, embed_coord(:,1), '.')

%%

figure;
plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.')
xlabel('rank from time')
ylabel('rank from vdm')


%%

figure;
scatter(embed_coord(:,1),embed_coord(:,2),50, time, '.')