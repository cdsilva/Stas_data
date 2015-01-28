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

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(images, ang_dis);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);

images_analyzed = order_all_images(images_registered, embed_coord);

plot_images(images_analyzed, dim)

%%
time = (1:nimages)';

corr(time, embed_coord,'type','spearman')

nsubimages = 10;
if corr(time, embed_coord) < 0
    embed_coord = -embed_coord;
end

fontsize = 8;

% idx_to_plot = randsample(nimages, nsubimages);
idx_to_plot = round(linspace(1, nimages, nsubimages));
[~, I] = sort(embed_coord(idx_to_plot));

make_fig(12, 12/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    imshow(make_gray_nuclei(images(:,:,:,idx_to_plot(j))));
end
saveas(gcf, 'drosophila_fixed_images_scrambled.pdf');


make_fig(12, 12/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    imshow(make_gray_nuclei(imrotate(images_registered(:,:,:,idx_to_plot(I(j))), -5, 'crop')));
    text(npixels/2, npixels/2, sprintf('%d', time(idx_to_plot(I(j)))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 4)
end
saveas(gcf, 'drosophila_fixed_images_ordered.pdf');

make_fig(12, 12/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    im_tmp = uint8(mean(images_analyzed(:,:,:,(j-1)*nimages/nsubimages+1:j*nimages/nsubimages), 4));
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




