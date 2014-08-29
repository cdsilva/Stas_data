clear all
close all

%%

npixels = 100;

dt = 0.5;

nmovies = 5;

nrot = 36;
nshifts = 0;
shift_max = 0.1;

neigs = 10;
alpha = 0;

channel = 1;

theta_err = zeros(nmovies, 1);
rank_corr = zeros(nmovies, 1);

%%

k = 1;

obj = VideoReader(sprintf('live_imaging/expt%d.mat', k));
image_set_full = obj.read;
nimages = get(obj, 'NumberOfFrames');
time = (1:nimages)' * dt;

image_set = zeros(npixels, npixels, nimages, 'uint8');
for i=1:nimages
    image_set(:, :, i) = image_fn(images(:, :, channel, i), npixels);
end

%
rng(21*k);
theta = 360*rand(nimages, 1);

for i=1:nimages
    image_set(:, :, i) = imrotate(image_set(:, :, i), theta(i), 'crop');
end

%
[R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
dim = size(R, 1) / nimages;

%
W2 = W.^2;
eps = median(W2(:))/10;
[R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);

%
theta_opt = zeros(size(theta));
for i=1:nimages
    R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
    theta_opt(i) = atan2d(R_tmp(2,1), R_tmp(1,1));
end
theta_err(k) =  min(std(mod(theta - theta_opt, 360)),std(mod(theta - theta_opt+180, 360)));

%
if corr(time, embed_coord(:,1)) < 0
    embed_coord(:,1) = -embed_coord(:,1);
end

rank_corr(k) = corr(compute_ranks(time), compute_ranks(embed_coord(:,1)));


[~, max_idx] = max(mod(theta, 360));
make_fig(4.3,4);
plot(mod(theta, 360), mod(theta_opt-theta_opt(max_idx)-10, 360)+10, '.k')
xlabel('true angle')
ylabel('recovered angle')

figure;
plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.k')
xlabel('true rank')
ylabel('recovered rank')

nprint_images = 10;
subplot_dim1 = nprint_images;
subplot_dim2 = 1;
plot_idx = round(linspace(3, nimages, nprint_images));

figure;
for j=1:nprint_images
    make_subplot(subplot_dim1, subplot_dim2, 0.01, j);
    imshow(image_set(:, :, plot_idx(j)));
end

figure;
[~, I] = sort(embed_coord(plot_idx, 1));
for j=1:nprint_images
    im_tmp = imrotate(image_set(:, :, plot_idx(I(j))), -theta_opt(plot_idx(I(j)))-75, 'crop');
    make_subplot(subplot_dim1, subplot_dim2, 0.01, j);
    imshow(im_tmp);
    text(npixels/2, npixels/2, sprintf('%2.1f min', time(plot_idx(I(j)))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
end

figure;
for j=1:nimages
    make_subplot(10, 4, 0.01, j);
    imshow(image_set(:, :, j));
end

figure;
[~, I] = sort(embed_coord(:, 1));
for j=1:nimages
    im_tmp = imrotate(image_set(:, :, I(j)), -theta_opt(I(j))-75, 'crop');
    make_subplot(10, 4, 0.01, j);
    imshow(im_tmp);
    text(npixels/2, npixels/2, sprintf('%2.1f min', time(I(j))),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
end



