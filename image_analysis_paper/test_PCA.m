clear all
close all

%res = '-r300';
%fmt = '-djpeg';

fontsize = 10;

im_save_dir = 'paper_figures2';

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
    %subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
    
end
saveas(gcf,sprintf('%s/raw_data1', im_save_dir), 'pdf')

image_set_aligned = zeros(size(image_set), 'uint8');
image_set_aligned_withnuclei = zeros(size(image_set), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    %         for j=1:3
    %             im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    %         end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set(:,:,:,i);
    for j=1:3
        im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned_withnuclei(:,:,:,i) = im1;
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

fontsize = 8;

data_reshaped = zeros(m, npixels*npixels*nchannels);

[eigenimages, D, proj_coeff] = PCA_images(image_set_aligned, 9);

figure;
for i=1:9
    subplot(3, 3, i)
    imshow(eigenimages(:, :, :, i))
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
imshow(eigenimages(:, :, :, 1),'border','tight')
saveas(gcf,sprintf('%s/PCA_eigenimage1', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
imshow(eigenimages(:, :, :, 2),'border','tight')
saveas(gcf,sprintf('%s/PCA_eigenimage2', im_save_dir), 'pdf')

[~, I] = sort(proj_coeff(:, 1));
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set_aligned_withnuclei(:,:,:,I(i));
    
    %subplot(subplot_dim1, subplot_dim2, i);
    %subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
    
end
saveas(gcf,sprintf('%s/PCA_ordered', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
plot(mem_lengths, proj_coeff(:,1)/1000,'.')
xlabel('membrane thickness', 'fontsize', fontsize)
ylabel('projection onto PC1', 'fontsize', fontsize)
saveas(gcf,sprintf('%s/PCA_corr', im_save_dir), 'pdf')

ranks_from_PCA = compute_ranks(proj_coeff(:,1));

corr(ranks_from_membranes, ranks_from_PCA)

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 4])
set(gcf, 'paperposition',[0 0 8 4])
plot(proj_coeff(:,1)/1000, proj_coeff(:,2)/1000,'.')
xlabel('projection onto PC1', 'fontsize', fontsize)
ylabel('projection onto PC2', 'fontsize', fontsize)
saveas(gcf,sprintf('%s/PCA_12', im_save_dir), 'pdf')

%%
nstages = 7;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
    make_subplot(nstages, 1, 0.01, i);
    stage_indices = I(max(1, round((i-1)*m/nstages)+1):min(m,round(i*m/nstages)));
    im1 = uint8(mean(double(image_set_aligned_withnuclei(:,:,:,stage_indices)), 4));
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
saveas(gcf,sprintf('%s/average_trajectory_PCA', im_save_dir), 'pdf')


%%
W2 = squareform(pdist(proj_coeff(:, 1:2))).^2;
eps2 = median(W2(:));

[V, D] = dmaps(W2, eps2, 10);

figure;
plot(V(:,2),mem_lengths,'.')
xlabel('\phi_2')
ylabel('membrane thickness')

[~, I] = sort(V(:,2));
if corr(V(:,2), mem_lengths) < 0
    V(:,2) = -V(:,2);
end
figure;
for i=1:m
    im1 = image_set_aligned(:,:,:,I(i));
    
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end

%%
nstages = 7;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
    make_subplot(nstages, 1, 0.01, i);
    stage_indices = I(max(1, round((i-1)*m/nstages)+1):min(m,round(i*m/nstages)));
    im1 = uint8(mean(double(image_set_aligned_withnuclei(:,:,:,stage_indices)), 4));
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
saveas(gcf,sprintf('%s/average_trajectory_DMAPS', im_save_dir), 'pdf')




