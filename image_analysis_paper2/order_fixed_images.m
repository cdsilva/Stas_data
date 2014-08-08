clear all
close all

%%

nimages = 132;
make_subplot = @(i) subplot(11, 12, i);

npixels = 100;
image_set = zeros(npixels, npixels, 3, nimages, 'uint8');
edge_tol = 1e3;
channel = 1;

for i=1:nimages
    image = imread(sprintf('../membrane_pictures/14_0501_dpERK_late/emb%02d.tif', i));
    image = imresize(image, [npixels npixels]);
    image(:, :, channel) = 0.5*adapthisteq(image(:, :, channel));
    image = crop_image(image, edge_tol, channel);
        
    image_set(:, :, :, i) = image;
end

%%

figure;
for i=1:nimages
    make_subplot(i)
    imshow(image_set(:, :, :, i))
end

%%

nrot = 20;
nshifts = 0;
shift_max = 0.1;
[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);

dim = size(R, 1) / nimages;

%%
W2 = W.^2;
eps = median(W2(:))/5;
neigs = 10;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs);

figure; plot(embed_coord(:,1), embed_coord(:,2),'.')


%%
figure;
for i=1:nimages
    make_subplot(i)
    imshow(rotate_image(image_set(:, :, :, i), R_opt(dim*(i-1)+1:dim*i, :)'))
end

%%
[~, I] = sort(embed_coord(:, 1));

figure;
for i=1:nimages
    make_subplot(i)
    imshow(rotate_image(image_set(:, :, :, I(i)), R_opt(dim*(I(i)-1)+1:dim*I(i), :)'))
end




