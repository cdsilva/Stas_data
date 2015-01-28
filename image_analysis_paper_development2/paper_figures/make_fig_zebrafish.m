clear all
close all

%%

% addpath ../code_for_distribution

%%
image_name = 'image';
image_ext = 'tif';
stack_name = '';
nimages = 120;
nstack = 0;
npixels = 100;
dim = 2;

channel_weight = 1;
channel_blur = 0;
channel_normalize = 0;
channel_mean_center = 1;
resize_image = false;

ang_dis = 10;

eps_scale = 0.25;

image_dir = '../code_for_distribution/zebrafish';

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

fid = fopen(sprintf('%s/times.txt', image_dir), 'r');
time = fscanf(fid, '%f');
fclose(fid);

figure;
plot(tiedrank(time), tiedrank(embed_coord),'.')

corr(time, embed_coord,'type','spearman')

fid = fopen(sprintf('%s/angles.txt', image_dir), 'r');
angles = fscanf(fid, '%f');
fclose(fid);

angles_opt = R_to_theta(R_opt);
[~, angles_opt] = angle_corr(angles, angles_opt);

%%

nsubimages = 10;
if corr(time, embed_coord) < 0
    embed_coord = -embed_coord;
end
% [~, I] = sort(embed_coord((1:nsubimages)+10));
% I = I+10;
[~, I] = sort(embed_coord);
I = I([1 20 35 50 60 70 80 94 100 120]);
I2 = sort(I);

make_fig(13, 13/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    imshow(images(:,:,I2(j)));
end
saveas(gcf, 'zebrafish_scrambled.pdf');

make_fig(13, 13/nsubimages);
for j=1:nsubimages
    make_subplot(nsubimages, 1, 0.01, j);
    imshow(imrotate(images_registered(:,:,I(j)), 70, 'crop'));
    text(npixels, npixels, sprintf('%2.0f min', time(I(j))-min(time)),'color',0.95*ones(1,3),'HorizontalAlignment','right','VerticalAlignment','bottom', 'fontsize', 4)

%     text(npixels/2, npixels/2, sprintf('%2.0f min', time(I(j))-min(time)),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
end
saveas(gcf, 'zebrafish_ordered.pdf');

fontsize = 8;
make_fig(3.25,3.25);
plot(tiedrank(time), tiedrank(embed_coord),'.k')
xlabel('true rank', 'fontsize', fontsize)
ylabel('recovered rank', 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gca, 'xtick', [0 50 100])
set(gca, 'ytick', [0 50 100])
axis([0 125 0 125])
axis square
saveas(gcf, 'zebrafish_rank_corr.pdf');

make_fig(4.3, 4);
plot(angles, angles_opt,'.k')
xlabel('true angle')
ylabel('recovered angle')
set(gca, 'xtick', [0 180 360]);
set(gca, 'ytick', [0 180 360]);
set(gca, 'xticklabel', {'0°'; '180°'; '360°'});
set(gca, 'yticklabel', {'0°'; '180°'; '360°'});
axis([0 400 0 400])
axis square

%% bootstrap
% nsamples = 10;
% nreps = 50;
% sample_vec = round(linspace(10, nimages, nsamples));
% theta_err = zeros(nsamples, nreps);
% rank_corr = zeros(nsamples, nreps);
%
% for j = 1:nsamples
%
%     j
%
%     for j2 = 1:nreps
%
%         sample_idx = randi(nimages, sample_vec(j), 1);
%
%         R_idx = reshape([2*sample_idx'-1; 2*sample_idx'], [], 1);
%
% %         [R, W] = compute_pairwise_alignments(images(:,:,sample_idx), nrot);
%
%         % compute optimal rotations + embedding coordinates using vector diffusion maps
%         [R_opt, embed_coord, D2] = vdm(R(R_idx, R_idx), W(sample_idx, sample_idx), eps_scale, ncomps);
%
%         % register images using optimal rotations
% %         images_registered = register_all_images(images(:,:,sample_idx), R_opt);
%
% %         images_analyzed = order_all_images(images_registered, embed_coord);
%
%
%         rank_corr(j, j2) = abs(corr(time(sample_idx), embed_coord,'type','spearman'));
%
%         angles_opt = R_to_theta(R_opt);
%         angles_opt = angles_opt + mean(atan2d(sind(angles(sample_idx)-angles_opt), cosd(angles(sample_idx)-angles_opt)));
%         angles_opt = mod(angles_opt, 360);
%         idx = find(angles(sample_idx) - angles_opt > 180);
%         angles_opt(idx) = angles_opt(idx) + 360;
%         idx = find(angles(sample_idx) - angles_opt < -180);
%         angles_opt(idx) = angles_opt(idx) - 360;
%
%         theta_err(j, j2) = corr(angles(sample_idx), angles_opt);
%     end
% end
%
% figure;
% plot(sample_vec, median(theta_err, 2), '.')
% set(gca, 'ylim', [0 1])
%
% figure;
% plot(sample_vec, median(rank_corr, 2), '.')
% set(gca, 'ylim', [0 1])

