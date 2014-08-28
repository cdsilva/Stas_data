clear all
close all

%%

npixels = 100;

nimages = 132;
image_dir = '../membrane_pictures/14_0501_dpERK_late';
image_name = 'emb';

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

ind = setdiff(1:nimages, [32 91 28 77 46 116 5 17 1 2 3 29]);

image_set = image_set(:, :, :, ind);
nimages = length(ind);

dim1 = 8;
dim2 = 15;

%%
% make_fig(17, dim1*(17/dim2));
% for i=1:nimages
%     make_subplot(dim2, dim1, 0.01, i);
%     imshow(make_gray_nuclei(image_set(:,:,:,i)))
% end

%%

nrot = 36;
nshifts = 0;
shift_max = 0.1;

[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages;

%
W2 = W.^2;
eps = median(W2(:))/10;
[R_opt, embed_coord_true, D2, D] = vdm(R, W2, eps, neigs, alpha);

%
if embed_coord_true(1,1) < 0
    embed_coord_true(:,1) = -embed_coord_true(:,1);
end

%%


nrot = 36;
nshifts = 0;
shift_max = 0.1;

[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages;

%
W2 = W.^2;
eps = median(W2(:))/10;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);

%
if embed_coord(1,1) < 0
    embed_coord(:,1) = -embed_coord(:,1);
end

figure;
plot(compute_ranks(embed_coord_true(:,1)), compute_ranks(embed_coord(:,1)),'.')

corr(compute_ranks(embed_coord_true(:,1)), compute_ranks(embed_coord(:,1)))