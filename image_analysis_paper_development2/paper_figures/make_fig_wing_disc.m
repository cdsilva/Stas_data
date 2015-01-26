
clear all
close all

%%
image_dir = '../code_for_distribution/wing_disc_stacks';
image_name = 'disc';
image_ext = 'tif';
stack_name = 'z';
nimages = 46;
nstack = 21;
npixels = 100;
dim = 3;

channel_weight = [1 1 1];
channel_blur = [0 0 0];
channel_normalize = [0 0 0];
channel_mean_center = [0 0 1];
resize_image = false;

nrot = 36;

eps_scale = 0.5;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(squeeze(max(images, [], ndims(images)-1)), nrot);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);
images_registered = imrotate(images_registered, 90);

plot_images(images_registered, dim)

W = squareform(pdist(reshape(double(images_registered), [], nimages)')).^2;
[embed_coord, D2] = dm(W, eps_scale, ncomps);

images_analyzed = order_all_images(images_registered, embed_coord);

%%

plot_images(immultiply(images_analyzed, 1.5), dim)
plot_images(immultiply(permute(images_analyzed, [2 4 3 1 5]), 1.5), dim)
plot_images(immultiply(permute(images_analyzed, [4 1 3 2 5]), 1.5), dim)

%%

plot_wing_disc_projections(images_registered(:,:,:,:,2), []);
annotation('arrow', [0.1 0.2], [0.1 0.1], 'color', 0.95*ones(1,3))
annotation('arrow', [0.1 0.1], [0.1 0.2], 'color', 0.95*ones(1,3))
annotation('textbox', [0.15 0.05 0.1 0.1], 'string','x','color', 0.95*ones(1,3), 'edgecolor','none', 'verticalalignment','middle','horizontalalignment','left')
annotation('textbox', [0.05 0.15 0.1 0.1], 'string','y','color', 0.95*ones(1,3), 'edgecolor','none', 'verticalalignment','bottom','horizontalalignment','center')

annotation('arrow', [0.85 0.9], [0.1 0.1], 'color', 0.95*ones(1,3))
annotation('arrow', [0.85 0.85], [0.1 0.15], 'color', 0.95*ones(1,3))
annotation('textbox', [0.87 0.05 0.1 0.1], 'string','z','color', 0.95*ones(1,3), 'edgecolor','none', 'verticalalignment','middle','horizontalalignment','left')
annotation('textbox', [0.82 0.11 0.1 0.1], 'string','y','color', 0.95*ones(1,3), 'edgecolor','none', 'verticalalignment','bottom','horizontalalignment','center')

annotation('arrow',  [0.1 0.1], [0.85 0.89],'color', 0.95*ones(1,3))
annotation('arrow',  [0.1 0.15],[0.85 0.85], 'color', 0.95*ones(1,3))
annotation('textbox', [0.05 0.84 0.1 0.1], 'string','z','color', 0.95*ones(1,3), 'edgecolor','none', 'verticalalignment','bottom','horizontalalignment','center')
annotation('textbox', [0.12 0.8 0.1 0.1], 'string','x','color', 0.95*ones(1,3), 'edgecolor','none', 'verticalalignment','middle','horizontalalignment','left')

saveas(gcf, 'wing_disc_example.pdf');

%%
fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

if corr(time, embed_coord) < 0
    embed_coord = -embed_coord;
end

nsubimages = 10;
[~, I] = sort(embed_coord(20+(1:nsubimages)));
I = I + 20;

plot_wing_disc_projections(images_registered(:,:,:,:,I), time(I)-min(time));
saveas(gcf, 'wing_disc_ordered.pdf');

return


% make_fig(17, 3*17/nsubimages);
% for j=1:nsubimages
%     make_subplot(nsubimages, 3, 0.01, j);
%     imshow(max(images_registered(:,:,:,:,I(j)),[], 4))
%     make_subplot(nsubimages, 3, 0.01, nsubimages+j);
%     imshow(max(permute(images_registered(:,:,:,:,I(j)), [1 4 3 2]),[], 4))
% 
%     make_subplot(nsubimages, 3, 0.01, 2*nsubimages+j);
%     imshow(max(permute(images_registered(:,:,:,:,I(j)), [4 2 3 1]),[], 4))
%     text(npixels/2, npixels/2, sprintf('~%2.0f min', time(I(j))-min(time)),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
% end

