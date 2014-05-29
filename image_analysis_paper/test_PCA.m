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
% saveas(gcf,sprintf('%s/raw_data1', im_save_dir), 'pdf')

image_set_aligned = zeros(size(image_set), 'uint8');
image_set_aligned_withnuclei = zeros(size(image_set), 'uint8');
theta_adjust = -58;
for i=1:m
    im_tmp = image_set(:,:,:,i);
    %         for j=1:3
    %             im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    %         end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set(:,:,:,i);
    for j=1:3
        im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
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

[eigenimages, D_PCA, proj_coeff] = PCA_images(image_set_aligned, m);

diag(D_PCA) / sum(diag(D_PCA))

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
% saveas(gcf,sprintf('%s/PCA_eigenimage1', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
imshow(eigenimages(:, :, :, 2),'border','tight')
% saveas(gcf,sprintf('%s/PCA_eigenimage2', im_save_dir), 'pdf')

[~, I] = sort(proj_coeff(:, 1));
if corr(I, mem_lengths) < 0 
    I = flipud(I);
end
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
% saveas(gcf,sprintf('%s/PCA_ordered', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
plot(mem_lengths, proj_coeff(:,1)/1000,'.')
xlabel('membrane thickness', 'fontsize', fontsize)
ylabel('projection onto PC1', 'fontsize', fontsize)
% saveas(gcf,sprintf('%s/PCA_corr', im_save_dir), 'pdf')

ranks_from_PCA = compute_ranks(proj_coeff(:,1));

corr(ranks_from_membranes, ranks_from_PCA)

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 4])
set(gcf, 'paperposition',[0 0 8 4])
plot(proj_coeff(:,1)/1000, proj_coeff(:,2)/1000,'.')
xlabel('projection onto PC1', 'fontsize', fontsize)
ylabel('projection onto PC2', 'fontsize', fontsize)
% saveas(gcf,sprintf('%s/PCA_12', im_save_dir), 'pdf')

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
% saveas(gcf,sprintf('%s/average_trajectory_PCA', im_save_dir), 'pdf')


%% DMAPS on PCA
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
% saveas(gcf,sprintf('%s/average_trajectory_DMAPS', im_save_dir), 'pdf')

%% VDM

eps = median(W(:))/10;
neigs = 3*m;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);


figure; 
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
plot(abs(diag(D)),'.')
xlabel('k')
ylabel('|\lambda_k|')
% saveas(gcf,sprintf('%s/data1_evals', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
nbins = 25;
hist(diag(D), nbins)
set(gca, 'xlim', [-1 1])
grid on
set(gca, 'ytick', 0:3:30)
set(gca, 'xtick', -1:0.2:1)
xlabel('\lambda')
ylabel('Pr(\lambda)')
% hold on
% [nelements, centers] = hist(diag(D(7:end, 7:end)), nbins);
% [alpha, rad] = fit_semicircle(centers, nelements);
% plot(centers, wigner_semicircle(centers, rad, alpha))
saveas(gcf,sprintf('%s/data1_evals_dist', im_save_dir), 'pdf')

return

figure;  
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
scatter(embed_idx(1,:), embed_idx(2, :), 500, var(embed_coord),'.')
colorbar
xlabel('k')
ylabel('l')
% saveas(gcf,sprintf('%s/data1_coord_var', im_save_dir), 'pdf')


theta_adjust = -55;

image_set_aligned = zeros(size(image_set), 'uint8');
image_set_aligned_withnuclei = zeros(size(image_set), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    %         for j=1:3
    %             im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    %         end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set(:,:,:,i);
    for j=1:3
        im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
    end
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    image_set_aligned_withnuclei(:,:,:,i) = im1;
end


%% order using vdm

idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);

[~, I] = sort(embed_coord(:, idx));
if corr(embed_coord(:, idx), mem_lengths) < 0
    I = flipud(I);
end
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    
    imshow(image_set_aligned_withnuclei(:,:,:,I(i)));
end
% saveas(gcf,sprintf('%s/VDM_data1_ordered', im_save_dir), 'pdf')

corr(embed_coord(:, idx), mem_lengths, 'type', 'spearman')

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
% saveas(gcf,sprintf('%s/average_trajectory_VDM', im_save_dir), 'pdf')

%%

idx1 = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);
%idx2 = find(embed_idx(1,:) == 7 & embed_idx(2,:) == 3);
idx2 = find(embed_idx(1,:) == 10 & embed_idx(2,:) == 1);

figure; plot(embed_coord(:, idx1), embed_coord(:, idx2), '.')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.')
hold on
dist_in_dmaps = squareform(pdist(embed_coord(:, [idx1 idx2])));
draw_tol = 0.01;
im_delta = 0.01;
hold on
for i=1:m
    if i == 1 || min(dist_in_dmaps(i, 1:i-1)) > draw_tol
        image('cdata', image_set_aligned(:, :, :, i), 'xdata', [embed_coord(i,idx1) embed_coord(i,idx1)+im_delta], 'ydata', [embed_coord(i,idx2) embed_coord(i,idx2)+im_delta])
    end
end
axis equal