clear all
close all

%%

npixels = 100;

nimages = 132;
image_dir = '../membrane_pictures/14_0501_dpERK_late';
image_name = 'emb';

nrot = 40;
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

%%

figure;
for i=1:nimages
    make_subplot(10, 14, 0.01, i);
    imshow(image_set(:,:,:,i))
end


%%
% ind = setdiff(1:nimages, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 32 102 101 100]);
%
% image_set = image_set(:, :, :, ind);
% nimages = length(ind);

ind = setdiff(1:nimages, [32 91 28 77 46 116]);

image_set = image_set(:, :, :, ind);
nimages = length(ind);

dim1 = 7;
dim2 = 18;

%%
make_fig(17, dim1*(17/dim2));
for i=1:nimages
    make_subplot(dim2, dim1, 0.01, i);
    imshow(make_gray_nuclei(image_set(:,:,:,i)))
end
saveas(gcf, 'fixed_images_unregistered_unordered.pdf');

%%

[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages;

%%
W2 = W.^2;
eps = median(W2(:))/10;
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

%%

rot_angle = -35;

make_fig(17, dim1*(17/dim2));
for i=1:nimages
    make_subplot(dim2, dim1, 0.01, i);
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,I(i)));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
end
saveas(gcf, 'fixed_images_registered_ordered.pdf');

%%
nstages = 14;
make_fig(17, 17/nstages);
for i=1:nstages
    make_subplot(nstages, 1, 0.01, i);
    idx = nimages/nstages*(i-1)+1:nimages/nstages*i;
    im_tmp = mean(double(image_set_aligned(:,:,:,I(idx))), 4);
    im_tmp = make_gray_nuclei(uint8(im_tmp));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
end
saveas(gcf, 'fixed_images_average_trajectory.pdf');

%%

figure;
for i=1:nimages
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,I(i)));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
    pause(0.1)
    clf
end
