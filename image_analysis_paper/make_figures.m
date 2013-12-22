clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';
print_figures = true;

dpERK_data = '../membrane_lengths/oct16.mat';
dpERK_image_dir = '../membrane_pictures/membrane2/dpERK_staining';
dpERK_membrane_dir = '../membrane_pictures/membrane2';

%% load data
load(dpERK_data);
dpERK = dpERK_raw;

[m, n] = size(dpERK);

% scramble data alignments
dpERK_unaligned = zeros(size(dpERK));
rng(12345);
rand_offsets = zeros(m,1);
for i=1:m
    rand_offsets(i) = randi(n);
    dpERK_unaligned(i,:) = circshift(dpERK(i,:),[0 rand_offsets(i)]);
end

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);
image_set_membrane = zeros(npixels, npixels, m);

image_channel = 2;
membrane_channel = 2;

% image indices to plot
im_save_idx = [1,9,17,25,30,34,38,43,49];

figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', dpERK_image_dir, i));

    % resize image
    im1 = imresize(im1, [npixels npixels]);

    % extract relevent color from image
    im1 = im1(:,:,image_channel);

    %store image
    im1 = double(im1);
    image_set(:, :, i) = im1;
    
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', dpERK_membrane_dir, i));

    % resize image
    im1 = imresize(im1, [npixels npixels]);

    % extract relevent color from image
    im1 = im1(:,:,membrane_channel);

    % adjust image
    im1 = imadjust(im1);
    
    %store image
    im1 = double(im1);
    image_set_membrane(:, :, i) = im1;
end

%% plot images

image_idx = 44;

im1 = imread(sprintf('%s/emb%02d.tif',dpERK_image_dir, image_idx));

im1_dpERK = im1;
im1_dpERK(:,:,1) = 0;
im1_dpERK(:,:,3) = 0;

im1_DI = im1;
im1_DI(:,:,1) = im1_DI(:,:,3);
im1_DI(:,:,2) = 0;

im1_membrane = imread(sprintf('%s/emb%02d.tif',dpERK_membrane_dir, image_idx));

figure;
imshow(im1_dpERK)
if print_figures
    print('drosophila_dpERK', fmt, res)
end

figure;
imshow(im1_DI)
if print_figures
    print('drosophila_DI', fmt, res)
end

figure;
imshow(im1_membrane)
if print_figures
    print('drosophila_membrane', fmt, res)
end

%% plot circle and line profiles

r = 1;
theta = linspace(0, 2*pi, n+1);
theta = theta(1:end-1);
theta = theta - pi/4;

x = r * cos(theta);
y = r * sin(theta);

figure;
scatter(x, y, 2000, dpERK(image_idx,:), '.')
hold on
delta=0.1;
plot([x(1)-delta x(1)+delta], [y(1)+delta y(1)-delta], 'linewidth',6, 'color','k')
text(x(1)+delta+0.05, y(1)-delta-0.05, 'Open circle')
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
set(gca, 'visible', 'off')
if print_figures
    print('circle_profile', fmt, res)
end

figure;
imagesc(repmat(dpERK(image_idx,:), 5, 1))
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
set(gca, 'visible', 'off')
if print_figures
    print('line_profile', fmt, res)
end


%% plot scrambled data

% plot data
figure;
imagesc(dpERK)
ylabel('data point (unordered)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_unordered',fmt,res)
end

%% plot data ordered by membrane length

% plot data
figure;
imagesc(dpERK_sort)
ylabel('data point (ordered using membrane thickness)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_ordered_membrane',fmt,res)
end

%% PCA

make_pca_figures;

%% DMAPS

make_dmaps_figures;

%% alignment

make_1d_alignment_figures;

%% raw dpERK images-- synchronization

make_2d_alignment_figures;

%% raw membrane pictures-- synchronization

make_2d_membrane_alignment_figures;

%% raw dpERK images

make_2d_scattering_figures;

%% membrane images

make_2d_membrane_scattering_figures;
    
