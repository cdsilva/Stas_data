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

nsubsamples = 15;
subsample_idx = round(linspace(1, nimages, nsubsamples));
make_fig(17, 17/nsubsamples);
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
make_fig(17, 17/nsubsamples);
for i=1:nsubsamples
    make_subplot(nsubsamples, 1, 0.01, i);
    im_tmp = make_gray_nuclei(image_set_aligned(:,:,:,subsample_idx(I2(i))));
    imshow(imrotate(im_tmp, rot_angle, 'crop'))
end
saveas(gcf, 'fixed_images_registered_ordered_subsampled.pdf');

%%
nstages = 15;
frame_points = linspace(1, nimages, nstages);
window_eps = 3^2;
window_tol = 0.01;

make_fig(17, 17/nstages);
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
