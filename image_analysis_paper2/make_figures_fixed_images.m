clear all
close all

%%

npixels = 100;

nimages = 132;
image_dir = '../membrane_pictures/14_0501_dpERK_late';
image_name = 'emb';

nrot = 36;
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

ind = setdiff(1:nimages, [32 91 28 77 46 116 5 17 1 2 3 29]);

image_set = image_set(:, :, :, ind);
nimages = length(ind);

dim1 = 8;
dim2 = 15;

%%
make_fig(17, dim1*(17/dim2));
for i=1:nimages
    make_subplot(dim2, dim1, 0.01, i);
    imshow(make_gray_nuclei(image_set(:,:,:,i)))
end
saveas(gcf, 'fixed_images_unregistered_unordered.pdf');

nsubsamples = 10;
subsample_idx = round(linspace(2, nimages, nsubsamples));
make_fig(12, 12/nsubsamples);
for i=1:nsubsamples
    make_subplot(nsubsamples, 1, 0.01, i);
    imshow(make_gray_nuclei(image_set(:,:,:,subsample_idx(i))))
end
saveas(gcf, 'fixed_images_unregistered_unordered_subsampled.pdf');

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

if embed_coord(1,1) < 0
    embed_coord(:,1) = -embed_coord(:,1);
end

[~, I] = sort(embed_coord(:,1));

%%

rot_angle = 85;

make_fig(17, dim1*(17/dim2));
for i=1:nimages
    make_subplot(dim2, dim1, 0.01, i);
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,I(i)));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
end
saveas(gcf, 'fixed_images_registered_ordered.pdf');

[~, I2] = sort(embed_coord(subsample_idx,1));
make_fig(12, 12/nsubsamples);
for i=1:nsubsamples
    make_subplot(nsubsamples, 1, 0.01, i);
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,subsample_idx(I2(i))));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
end
saveas(gcf, 'fixed_images_registered_ordered_subsampled.pdf');

%%
nstages = 10;
frame_points = linspace(1, nimages, nstages);
window_eps = 3^2;
window_tol = 0.01;

make_fig(12, 12/nstages);
for i=1:nstages
    
    window_weights = exp(-(frame_points(i)-(1:nimages)).^2/window_eps);
    window_weights = window_weights / sum(window_weights);

    im_tmp = immultiply(image_set_aligned(:,:,:,I(1)), window_weights(1));
    for j=2:nimages
            im_tmp = imlincomb(1, im_tmp, window_weights(j), image_set_aligned(:,:,:,I(j)));
    end
    im_tmp = make_gray_nuclei(im_tmp);
    
    make_subplot(nstages, 1, 0.01, i);
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
    
end
saveas(gcf, 'fixed_images_average_trajectory.pdf');

%%

I_manual = [8 48 94    81 65  118  55    14 27   67    62    82    16 53 109     5      115        45        21       104        36 ...
    33    40  91 100   26 15    51    85    97  111  28    32    61   101   105       108   112    86    23    72    93   119 ...
    54    73     3    39    58    71    75    60       10    49    46    41    22    12    38    43    56    66    68    87    17    13 ...
    117   113    83    90    76    34   103   107    88    98    30    25    57    70    29    19    31    89    80 ...
   106    20    50    99     2        77    44    78    18       110  95    69     1    84   120     9    47 ...
     4    42    96    74    92    37   116       59     7   114    24 6 35   79    52    64   102    11    63];

make_fig(17, dim1*(17/dim2));
for i=1:nimages
    make_subplot(dim2, dim1, 0.01, i);
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,I_manual(i)));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
    title(sprintf('%d', I_manual(i)))
end

r1(I_manual) = 1:nimages;
r2(I) = 1:nimages;
make_fig(4.5,4.5);
plot(r1, r2,'.')
xlabel('manual rank')
ylabel('algorithm rank')
set(gca, 'xtick', [0 50 100])
set(gca, 'ytick', [0 50 100])
axis([0 130 0 130])
axis square
saveas(gcf, 'rank_corr_fixed_images.pdf');

return
%%

writerObj = VideoWriter('gastrulation.avi');
writerObj.FrameRate = 50;
open(writerObj);

nframes = 1000;
frame_points = linspace(1, nimages, nframes);
window_eps = 3^2;

make_fig(8,8);

i = 1;
window_weights = exp(-(frame_points(i)-(1:nimages)).^2/window_eps);
window_weights = window_weights / sum(window_weights);

im_tmp = immultiply(image_set_aligned(:,:,:,I(1)), window_weights(1));
for j=2:nimages
    im_tmp = imlincomb(1, im_tmp, window_weights(j), image_set_aligned(:,:,:,I(j)));
end
im_tmp = make_gray_nuclei(im_tmp);

imshow(imrotate(im_tmp, rot_angle, 'crop'))
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

clf

for i=1:nframes
    
    window_weights = exp(-(frame_points(i)-(1:nimages)).^2/window_eps);
    window_weights = window_weights / sum(window_weights);
    
    im_tmp = immultiply(image_set_aligned(:,:,:,I(1)), window_weights(1));
    for j=2:nimages
        im_tmp = imlincomb(1, im_tmp, window_weights(j), image_set_aligned(:,:,:,I(j)));
    end
    im_tmp = make_gray_nuclei(im_tmp);
    
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
    frame = getframe;
    writeVideo(writerObj,frame);
    clf
end

close(writerObj);
