clear all
close all

%%

npixels = 100;

nimages = 132;
image_dir = '../membrane_pictures/14_0501_dpERK_late';
image_name = 'emb';

nrot = 20;
nshifts = 0;
shift_max = 0.1;

neigs = 10;
alpha = 0;

channel = 1;
nchannels = 3;

image_set = zeros(npixels, npixels, nchannels, nimages, 'uint8');

%%

for i=1:nimages
    im_tmp = imread(sprintf('%s/%s%02d.tif', image_dir, image_name, i));
    image_set(:, :, :, i) = image_fn_color(im_tmp, channel, npixels);
end

% image_set(:,:,2:3,:) = 0;


%%

figure;
for i=1:nimages
    subplot(12,11,i)
    imshow(image_set(:,:,:,i))
end


%%
ind = setdiff(1:nimages, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 32 102 101 100]);

image_set = image_set(:, :, :, ind);
nimages = length(ind);

%%

[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages;

%
W2 = W.^2;
eps = median(W2(:))/20;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);

%%
image_set_aligned = zeros(npixels, npixels, nchannels, nimages, 'uint8');
for i=1:nimages
    R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
    image_set_aligned(:, :, :, i) = rotate_image(image_set(:,:, :,i), R_tmp');
end

%%

if embed_coord(1,1) > 0
    embed_coord(:,1) = -embed_coord(:,1);
end

[~, I] = sort(embed_coord(:,1));
figure;
for i=1:nimages
    subplot(6,22,i)
    imshow(image_set_aligned(:,:,:,I(i)))
end
