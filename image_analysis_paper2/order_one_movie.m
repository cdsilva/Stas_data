clear all
close all

%%

npixels = 100;
[images, time] = read_video('../live_imaging/bomyi_emb01_gast01.avi', npixels);

nimages = length(time);
make_subplot = @(i) subplot(10, 10, i);

image_set = zeros(npixels, npixels, nimages, 'uint8');
edge_tol = 1e3;
channel = 1;

for i=1:nimages
    image = images(:, :, channel, i);
    image = imresize(image, [npixels npixels]);
    image = adapthisteq(image);
    image = crop_image(image, edge_tol, channel);
        
    image_set(:, :, i) = image;
end

%%

figure;
for i=1:nimages
    make_subplot(i)
    imshow(image_set(:, :, i))
end

%%

rng(123);
theta = 360*rand(nimages, 1);
idx = randperm(nimages);

for i=1:nimages
    image_set(:, :, i) = imrotate(image_set(:, :, i), theta(i), 'crop');
end

image_set = image_set(:, :, idx);
time = time(idx);

figure;
for i=1:nimages
    make_subplot(i)
    imshow(image_set(:, :, i))
end

%%

nrot = 20;
nshifts = 0;
shift_max = 0.1;
[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);

dim = size(R, 1) / nimages;

%%
W2 = W.^2;
eps = median(W2(:));
neigs = 10;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs);

%%
figure;
for i=1:nimages
    make_subplot(i)
    imshow(rotate_image(image_set(:, :, i), R_opt(dim*(i-1)+1:dim*i, :)'))
end

%%
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
scatter(embed_coord(:,1),embed_coord(:,2),50, time, '.')