clear all
close all

%%

npixels = 100;

nmovies = 6;

idx_start = [21 31 27 24 18 17];
% idx_start = ones(1, 6);
idx_end = [56 60 47 44 46 42];

shift_times = [1.6302 2.8781 1.0172 -1.7739 -0.9879 -4.9915];
dt = 0.5;

images_all = [];
time_all = [];

for k=1:nmovies
    load(sprintf('movie%d.mat', k));
    idx = idx_start(k):idx_end(k);
    
    images_all = cat(4, images_all, images(:, :, :, idx));
    time_all = [time_all; dt*time(idx)-shift_times(k)];
end

time = time_all;
images = images_all;

%%
nimages = length(time);
make_subplot = @(i) subplot(ceil(nimages/10), 10, i);

image_set = zeros(npixels, npixels, nimages, 'uint8');
channel = 1;

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
% 
% figure;
% for i=1:nimages
%     make_subplot(i)
%     imshow(image_set(:, :, i))
% end

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

% figure;
% for i=1:nimages
%     make_subplot(i)
%     imshow(image_set(:, :, i))
% end

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
% figure;
for i=1:nimages
    R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
    theta_opt(i) = atan2d(R_tmp(2,1), R_tmp(1,1));

%     make_subplot(i)
%     imshow(rotate_image(image_set(:, :, i), R_opt(dim*(i-1)+1:dim*i, :)'))
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

% figure;
% scatter(embed_coord(:,1),embed_coord(:,2),50, time, '.')

