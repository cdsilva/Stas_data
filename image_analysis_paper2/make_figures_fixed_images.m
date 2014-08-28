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

images_to_remove = [32 91 28 77 46 116 5 17 1 2 3 29];
ind = setdiff(1:nimages, images_to_remove);

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

% I_manual = [17 29 92 64 5 76 71 35 12 19 74 57 93 21 105  121 54 108 130 2 ...
% 1 62 3 44 127 20 36 9 111 49 96 120 41 26 131 112 60 123 30 70 ...
% 102 63 34 117 124 115 104 82 48 97 83 40 81 69 85 67 7 14 58 50 ...
% 55 27 65 22 129 75 78 18 52 119 16 94 98 125 47 109 101 86 37 114 ...
% 33 42 110 79 59 38 66 99 24 4 6 39 13 87 90 88 8 106 103 132 ...
% 56 95 118 80 23 100 68 45 107 126 31 25 128 51 84 53 11 43 10 73 ...
% 77 32 116 46 122 113 61 89 91 15 72 28];

I_manual = [17 5 29 64 92 35 12 71 57 76 19 74 93 21 105 54 121 130 108 127 ...
2 3 9 44 62 1 26 20 49 36 111 131 120 123 112 117 96 60 70 30 ...
41 124 115 102 34 104 63 48 97 82 81 83 40 69 85 67 58 7 14 50 ...
55 27 65 129 18 22 75 52 78 16 119 47 98 94 125 101 86 114 33 37 ...
109 42 100 110 80 99 38 59 66 24 6 4 39 90 56 8 13 106 88 87 ...
132 95 103 118 68 107 31 25 23 79 45 126 53 11 128 51 84 10 43 73 ...
77 32 116 46 122 15 113 61 89 91 72 28];

for i=sort(images_to_remove, 'descend');
    I_manual(I_manual == i) = [];
    idx = find(I_manual > i);
    I_manual(idx) = I_manual(idx) - 1;
end

make_fig(17, dim1*(17/dim2));
for i=1:nimages
    make_subplot(dim2, dim1, 0.01, i);
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,I_manual(i)));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
    title(sprintf('%d', I_manual(i)))
end

rank_manual(I_manual) = 1:nimages;
rank_algorithm(I) = 1:nimages;
make_fig(4.5,4.5);
plot(rank_manual, rank_algorithm,'.')
xlabel('expert rank')
ylabel('algorithm rank')
set(gca, 'xtick', [0 50 100])
set(gca, 'ytick', [0 50 100])
axis([0 130 0 130])
axis square
saveas(gcf, 'rank_corr_fixed_images.pdf');

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
    text(npixels/2, npixels/2, sprintf('%d', rank_manual(subsample_idx(I2(i)))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
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
