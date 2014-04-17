clear all
close all

%res = '-r300';
%fmt = '-djpeg';

fontsize = 10;

%im_save_dir = 'paper_figures';

%% directories where things are stored
time_data = '../membrane_lengths/oct16.mat';
image_dir = '../membrane_pictures/membrane2/dpERK_staining';

%% load membrane lengths

load(time_data);
mem_lengths = L(:,1);

% remove some embryos
ind = setdiff(1:52, [9, 15, 22]);

m = length(ind);

mem_lengths = mem_lengths(ind);

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% image parameters

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim1:-1:1)/subplot_dim1, (1:subplot_dim2)/subplot_dim2);

nchannels = 3;
image_set = zeros(npixels, npixels, nchannels, m, 'uint8');
image_set_raw = zeros(1024, 1024, nchannels, m, 'uint8');

nuclei = zeros(npixels, npixels, m, 'uint8');

%% load images

for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    if i==1
        im1 = circshift(im1, [0 150 0]);
    end
    image_set_raw(:, :, :, i) = im1;
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    nuclei(:, :, i) = im1(:, :, 1);
    
    im1(:, :, 1) = im1(:, :, 2);
    im1(:, :, 2) = im1(:, :, 3);
    im1(:, :, 3) = 0;
    
    %store image
    image_set(:, :, :, i) = im1;
    
end

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 10;
shift_step = 2;
dim = 3;

%[R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);

%save('pairwise_alignments_small.mat', 'R', 'W');
%load('pairwise_alignments_small.mat');

%save('pairwise_alignments_small_nopbc.mat', 'R', 'W');
load('pairwise_alignments_small_nopbc.mat');


%% synchronization

R_opt = ang_synch(R, dim);

[~, D] = eig(R);

figure;
hist(diag(D), 200)

%%

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set(:,:,:,i);
    
    for j=1:3
        im1(:,:,j) = im1(:,:,j) + nuclei(:,:,i);
    end
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end

image_set_aligned = zeros(size(image_set), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    %     for j=1:3
    %         im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    %     end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set_aligned(:,:,:,i);
    
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end

%%

npixels2 = 20;

data_reshaped = zeros(m, npixels2*npixels2*nchannels);

figure;
for i=1:m
    im1 = image_set_aligned(:,:,:,i);
    im1 = imresize(im1, [npixels2 npixels2]);
    data_reshaped(i, :) = im1(:)';
    
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end

data_mean = mean(data_reshaped);

data_reshaped = data_reshaped - repmat(data_mean, m, 1);

[V, D] = PCA(data_reshaped, 9);

figure;
bar(diag(D))
xlabel('k')
ylabel('\lambda_k (PCA)')

figure;
im1 = V(:,1);
im1 = im1 - min(im1);
im1 = im1 * 255 / max(im1);
im1 = reshape(im1, npixels2, npixels2, nchannels);
imshow(uint8(im1));

figure;
plot(data_reshaped*V(:,1),mem_lengths,'.')
xlabel('\langle x, V_1 \rangle')
ylabel('membrane thickness')


figure;
plot(data_reshaped*V(:,1),data_reshaped*V(:,2),'.')
xlabel('\langle x, V_1 \rangle')
ylabel('\langle x, V_2 \rangle')

