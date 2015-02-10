
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

ang_dis = 10;

eps_scale = 0.5;

[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, dim);
plot_images(images_raw, dim)

images = apply_image_functions(images_raw, npixels, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image);
plot_images(images, dim)

[R, W] = compute_pairwise_alignments(squeeze(max(images, [], ndims(images)-1)), ang_dis);

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
% [R_opt, embed_coord, D2] = vdm(R, ones(nimages), eps_scale, ncomps);
R_opt = ang_synch(R, 2);

% register images using optimal rotations
images_registered = register_all_images(images, R_opt);
images_registered = imrotate(images_registered, 90);
images_registered = images_registered(:,end:-1:1,:,:,:);

plot_images(images_registered, dim)

W = squareform(pdist(reshape(double(images_registered), [], nimages)')).^2;
[embed_coord, D2] = dm(W, eps_scale, ncomps);

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

if corr(time, embed_coord) < 0
    embed_coord = -embed_coord;
end

images_analyzed = order_all_images(images_registered, embed_coord);

%%

plot_images(immultiply(images_analyzed, 1.5), dim)
plot_images(immultiply(permute(images_analyzed, [2 4 3 1 5]), 1.5), dim)
plot_images(immultiply(permute(images_analyzed, [4 1 3 2 5]), 1.5), dim)

%%

fontsize = 6;

make_fig(6, 6)
% set(gcf, 'Renderer', 'Painters')

plot_wing_disc_projections(images_registered(:,:,:,:,2), [], 1, 1)
annotation('arrow', [0.05 0.15], [0.05 0.05], 'color', 0.95*ones(1,3), 'headlength', 5, 'headwidth', 5)
annotation('arrow', [0.05 0.05], [0.05 0.15], 'color', 0.95*ones(1,3), 'headlength', 5, 'headwidth', 5)
annotation('textarrow', [0.15 0.05], [0.05 0.05], 'color', 0.95*ones(1,3), 'string','x', 'headstyle','none', 'fontsize', fontsize)
annotation('textarrow', [0.05 0.05], [0.15 0.05], 'color', 0.95*ones(1,3), 'string','y', 'headstyle','none', 'fontsize', fontsize)

annotation('arrow', [0.82 0.87], [0.05 0.05], 'color', 0.95*ones(1,3), 'headlength', 5, 'headwidth', 5)
annotation('arrow', [0.82 0.82], [0.05 0.1], 'color', 0.95*ones(1,3), 'headlength', 5, 'headwidth', 5)
annotation('textarrow', [0.87 0.82], [0.05 0.05], 'color', 0.95*ones(1,3), 'string','z', 'headstyle','none', 'fontsize', fontsize)
annotation('textarrow', [0.82 0.82], [0.1 0.05], 'color', 0.95*ones(1,3), 'string','y', 'headstyle','none', 'fontsize', fontsize)

annotation('arrow',  [0.05 0.05], [0.82 0.87],'color', 0.95*ones(1,3), 'headlength', 5, 'headwidth', 5)
annotation('arrow',  [0.05 0.1],[0.82 0.82], 'color', 0.95*ones(1,3), 'headlength', 5, 'headwidth', 5)
annotation('textarrow',  [0.05 0.05], [0.87 0.82],'color', 0.95*ones(1,3), 'string','z', 'headstyle','none', 'fontsize', fontsize)
annotation('textarrow',  [0.1 0.05],[0.82 0.82], 'color', 0.95*ones(1,3), 'string','x', 'headstyle','none', 'fontsize', fontsize)

annotation('doublearrow',  [0 0.75], [0.76 0.76],'color', 'k', 'headlength', 5, 'headwidth', 5)
annotation('textbox', [0.75/2-0.25 0.77 0.5 0.2], 'string', 'A-P axis; 380 \mum', 'fontsize', 5, 'VerticalAlignment','bottom','HorizontalAlignment','center', 'Margin',0, 'EdgeColor','none')

annotation('doublearrow',  [0.8 0.9], [0.76 0.76],'color', 'k', 'headlength', 5, 'headwidth', 5)
annotation('textbox', [0.85-0.25 0.77 0.5 0.2], 'string', '21 \mum', 'fontsize', 5, 'VerticalAlignment','bottom','HorizontalAlignment','center', 'Margin',0, 'EdgeColor','none')

annotation('doublearrow', [0.76 0.76], [0.8 0.9],'color', 'k', 'headlength', 5, 'headwidth', 5)
annotation('textarrow', [0.79 0.79], [0.85 0.85], 'string', '21 \mum', 'fontsize', 5, 'VerticalAlignment','bottom','HorizontalAlignment','center', 'textrotation',-90, 'linestyle','none','headstyle','none','color', 0.99*ones(1,3), 'textcolor','k')

annotation('doublearrow', [0.76 0.76], [0 0.75],'color', 'k', 'headlength', 5, 'headwidth', 5)
annotation('textarrow', [0.79 0.79], [0.75/2 0.75/2], 'string', 'D-V axis; 380 \mum', 'fontsize', 5, 'VerticalAlignment','bottom','HorizontalAlignment','center', 'textrotation',-90, 'linestyle','none','headstyle','none','color', 0.99*ones(1,3), 'textcolor','k')

print('wing_disc_example.eps', '-depsc', '-r300');

%%

nsubimages = 15;
[~, I] = sort(embed_coord(20+(1:nsubimages)));
I = I + 20;

make_fig(10, 6);
plot_wing_disc_projections(images_registered(:,:,:,:,I), time(I), 5, 3)
% plot_wing_disc_projections(images_registered(:,:,:,:,I), time(I)-min(time));
% saveas(gcf, 'wing_disc_ordered.pdf');
print('wing_disc_ordered.eps', '-depsc', '-r300');

%%

nsubimages = 6;

avg_images = compute_average_trajectory(images_analyzed, nsubimages, 4);
make_fig(17, 17/nsubimages);
plot_wing_disc_projections(avg_images, [], nsubimages, 1)
print('wing_disc_average.eps', '-depsc', '-r300');

%%

ncomps = 20;
[embed_coord, D2] = dm(W, eps_scale, ncomps);

make_fig(8,8)
plot(abs(D2), '.')
xlabel('embedding coordinate')
ylabel('product of eigenvalues')
axis square
saveas(gcf, 'wing_disc_eval_spectrum.pdf');


