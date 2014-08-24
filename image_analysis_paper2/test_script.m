clear all
close all

%%

image = imread('../membrane_pictures/14_0501_dpERK_late/emb01.tif');

nimages = 20;
npixels = 100;

image = imresize(image, [npixels npixels]);

image_set = zeros(npixels, npixels, 3, nimages, 'uint8');

for i=1:nimages
    im_tmp = imrotate(image, 360*rand, 'crop');
%     im_tmp = circshift(im_tmp, randi([-5 5], 1, 2));
    image_set(:,:,:,i) = im_tmp;
end

%%

figure;
for i=1:nimages
    subplot(4,5,i)
    imshow(image_set(:,:,:,i))
end

%%

nrot = 10;
nshifts = 0;
shift_max = 0.05;

[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages;

%%
W2 = W.^2;
eps = median(W2(:))/10;
neigs = 10;
alpha = 0;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);

%%
% image_set_aligned = zeros(size(image_set), 'uint8');
% for i=1:nimages
%     R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
%     image_set_aligned(:, :, :, i) = rotate_image(image_set(:,:, :,i), R_tmp');
% end

image_set_aligned = register_all_images(image_set, R_opt);

%%

figure;
for i=1:nimages
    subplot(4,5,i)
    imshow(image_set_aligned(:,:,:,i))
end