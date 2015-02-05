clear all
close all

%%

% addpath ../code_for_distribution

%%
image_name = 'emb';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
npixels = 100;
dim = 2;

channel_weight = [0.5 1 1];
channel_blur = [0.05 0.05 0.05];
channel_normalize = [1 0 0];
channel_mean_center = [1 0 0];
resize_image = true;

ang_dis = 10;

eps_scale = 0.25;

image_dir = '../code_for_distribution/drosophila_fixed_images';

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, npixels, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(images, ang_dis);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

time = (1:nimages)';
if corr(time, embed_coord) < 0
    embed_coord = -embed_coord;
end
corr(time, embed_coord,'type','spearman')

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

images_analyzed = order_all_images(images_registered, embed_coord);

plot_images(images_analyzed, dim)

%%

fontsize = 8;
nsubimages = 10;

% idx_to_plot = randsample(nimages, nsubimages);
idx_to_plot = round(linspace(1, nimages, nsubimages));
[~, I] = sort(embed_coord(idx_to_plot));

make_fig(12, 12/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    imshow(make_gray_nuclei(images(:,:,:,idx_to_plot(j))));
    if j == 1
        hold on
        plot([0.1*npixels (0.1+0.8/3)*npixels], [0.05*npixels 0.05*npixels], 'color', 0.95*ones(1,3))
    end
end
saveas(gcf, 'drosophila_fixed_images_scrambled.pdf');

make_fig(12, 12/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    imshow(make_gray_nuclei(imrotate(images_registered(:,:,:,idx_to_plot(I(j))), -5, 'crop')));
    text(npixels/2, npixels/2, sprintf('%d', time(idx_to_plot(I(j)))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 8)
end
saveas(gcf, 'drosophila_fixed_images_ordered.pdf');

make_fig(12, 12/nsubimages);
avg_images = compute_average_trajectory(images_analyzed, nsubimages, 4);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
%     im_tmp = uint8(mean(images_analyzed(:,:,:,(j-1)*nimages/nsubimages+1:j*nimages/nsubimages), 4));
    im_tmp = avg_images(:,:,:,j);
    imshow(make_gray_nuclei(imrotate(im_tmp, -5, 'crop')));
end
saveas(gcf, 'drosophila_fixed_images_average.pdf');


make_fig(4,4);
plot(tiedrank(time), tiedrank(embed_coord),'.k')
xlabel('expert rank', 'fontsize', fontsize)
ylabel('algorithm rank', 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gca, 'xtick', [0 50 100])
set(gca, 'ytick', [0 50 100])
axis([0 125 0 125])
axis square
saveas(gcf, 'drosophila_fixed_images_rank_corr.pdf');

%%

ncomps = 20;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

make_fig(8,8)
plot(abs(D2), '.')
xlabel('embedding coordinate')
ylabel('product of eigenvalues')
axis square
saveas(gcf, 'drosophila_fixed_eval_spectrum.pdf');

%%

[Y, X] = meshgrid(1:npixels, 1:npixels);
X = X - (npixels+1)/2;
Y = Y - (npixels+1)/2;

theta = atan2d(Y, X);

nslices = 36;
theta_vec = linspace(-180, 180, nslices + 1);
twi_peaks = zeros(nimages, nslices);
for i=1:nslices
    idx = find(theta < theta_vec(i+1) & theta >= theta_vec(i));
    for j=1:nimages
        im_tmp = images_registered(:,:,3,j);
        twi_peaks(j, i) = mean(im_tmp(idx));
    end
end

std((twi_peaks * theta_vec(1:end-1)')./sum(twi_peaks, 2))













