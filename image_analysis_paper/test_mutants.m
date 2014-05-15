clear all
close all

%% directories where things are stored

image_dir = '../membrane_pictures/snanull_big';
m = 98;

% image_dir = '../membrane_pictures/snanull_old';
% m = 54;

%image_dir = '../membrane_pictures/snanull_old2';
%m = 47;

ind = 1:m;

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

figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    image_set_raw(:, :, :, i) = im1;
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    nuclei(:, :, i) = im1(:, :, 1);
    
    im1(:, :, 1) = im1(:, :, 2);
    im1(:, :, 2) = im1(:, :, 3);
    im1(:, :, 3) = 0;
    
    %store image
    image_set(:, :, :, i) = im1;
    
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
end

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 10;
shift_step = 2;
dim = 3;

%matlabpool open 2;
%[R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
%matlabpool close

%save('pairwise_alignments_snanull.mat', 'R', 'W');
load('pairwise_alignments_snanull.mat');


%% select which images to use

% ind = setdiff(1:m, [61 93 39 9 76 89 40 63 83 58]);
%
% image_set = image_set(:, :, :, ind);
% image_set_raw = image_set_raw(:, :, :, ind);
% nuclei = nuclei(:, :, ind);
%
% W = W(ind, ind);
%
% R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
% R = R(R_ind, R_ind);
%
% m = length(ind);
%
% subplot_dim1 = floor(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);

%% just dmaps

% for i=1:m
%     W(i, i) = 0;
% end

eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

figure;
plot(V(:,2),V(:,3),'.')


%% synchronization

R_opt = ang_synch(R, dim);

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

%% dmaps after synch

figure;
for i=1:m
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    
    imshow(image_set_aligned(:,:,:,i));
end

W2 = inf(m);

for i=1:m
    im1 = image_set_aligned(:,:,:,i);
    for j=1:i-1
        im2 = image_set_aligned(:,:,:,j);
        W2(i, j) = sum((double(im1(:)) - double(im2(:))).^2);
        W2(j, i) = W2(i, j);
    end
end

eps2 = median(W2(:));

[V2, D] = dmaps(W2, eps2, 10);

figure;
plot(V2(:,2),V2(:,3),'.')


figure;
im_delta = 0.02;
plot(V2(:,2),V2(:,3),'.')
hold on
for i=1:m
    image('cdata', image_set_aligned(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
    hold on
end

return

%% VDM

eps = median(W(:));
neigs = 10;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

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
%if corr(embed_coord(:, idx), mem_lengths) < 0
%    I = flipud(I);
%end
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    
    imshow(image_set_aligned_withnuclei(:,:,:,I(i)));
end

%%
nstages = 14;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    make_subplot(nstages, 1, 0.01, i);
    stage_indices = I(max(1, round((i-1)*m/nstages)+1):min(m,round(i*m/nstages)));
    im1 = uint8(mean(double(image_set_aligned_withnuclei(:,:,:,stage_indices)), 4));
    
    imshow(im1,'initialmagnification','fit','border','tight')
end

%%

eps = median(W(:));
[V2, D2] = dmaps(W, eps, 10);

figure;
plot(V2(:,2),V2(:,3),'.')

figure;
im_delta = 0.02;
plot(V2(:,2),V2(:,3),'.')
hold on
for i=1:m
    image('cdata', image_set_aligned_withnuclei(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
    hold on
end


%%

embed_coord2 = embed_coord;
for i=1:size(embed_coord, 2)
    for j=1:2
        embed_coord2(:, i) = embed_coord2(:, i) * sqrt(D(embed_idx(j,i), embed_idx(j,i)));
    end
end
W2 = squareform(pdist(embed_coord2)).^2;
eps2 = median(W2(:))*100;

[V2, D2] = dmaps(W2, eps2, 10);

figure;
plot(V2(:,2),V2(:,3),'.')




