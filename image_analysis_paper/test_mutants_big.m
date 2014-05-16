clear all
close all

%% directories where things are stored

im_save_dir = 'paper_figures2';

image_dir = '../membrane_pictures/snanull_big';
m = 98;

wt = [1 2 3 4 5 6 8 9 10 12 14 19 22 23 24 25 26 27 28 30 ...
    31 35 39 42 43 44 45 46 47 50 51 58 59 64 65 66 67 69 70 75 ...
    76 77 78 79 80 88 89 91 92 94 95]; %51

mut = [7 11 13 15 16 17 18 20 21 29 32 33 34 36 37 38 40 41 48 49 ...
    52 53 54 55 56 57 60 61 62 63 68 71 72 73 74 81 82 83 84 85 ...
    86 87 90 93 96 97 98]; %47

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
    im1(:, :, 2) = uint8(0.5*double(im1(:, :, 3)));
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

matlabpool open 4;
[R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
matlabpool close

%save('pairwise_alignments_snanull.mat', 'R', 'W');
% load('pairwise_alignments_snanull.mat');


%% select which images to use

ind = setdiff(1:m, [17]);

wt = [1 2 3 4 5 6 8 9 10 12 14 19 22 23 24 25 26 27 28 30 ...
    31 35 39 42 43 44 45 46 47 50 51 58 59 64 65 66 67 69 70 75 ...
    76 77 78 79 80 88 89 91 92 94 95]; %51

mut = [7 11 13 15 16 18 20 21 29 32 33 34 36 37 38 40 41 48 49 ...
    52 53 54 55 56 57 60 61 62 63 68 71 72 73 74 81 82 83 84 85 ...
    86 87 90 93 96 97 98]; %47

i = find(wt > 17);
wt(i) = wt(i) - 1;
i = find(mut > 17);
mut(i) = mut(i) - 1;


image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);
nuclei = nuclei(:, :, ind);

W = W(ind, ind);

R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
R = R(R_ind, R_ind);

m = length(ind);

subplot_dim1 = floor(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

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

W2 = zeros(m);
for i=1:m
    im1 = image_set_aligned(:,:,:,i);
    for j=1:i-1
        im2 = image_set_aligned(:,:,:,j);
        W2(i, j) = sum((double(im1(:)) - double(im2(:))).^2);
        W2(j, i) = W2(i, j);
    end
end

eps2 = median(W2(:))*10;

[V2, D] = dmaps(W2, eps2, 10);

figure;
plot(V2(:,2),V2(:,3),'.')


figure;
im_delta = 0.04;
plot(V2(:,2),V2(:,3),'.')
hold on
for i=1:m
    image('cdata', image_set_aligned(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
    hold on
end
axis equal

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
im_delta = 0.04;
plot(V2(:,2),V2(:,3),'.')
dist_in_dmaps = squareform(pdist(V2(:, 2:3)));
draw_tol = 0.04;
hold on
for i=1:m
    if i == 1 || min(dist_in_dmaps(i, 1:i-1)) > draw_tol
        image('cdata', image_set_aligned(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
        %hold on
    end
end
axis equal


%%
ranks_from_dmaps = compute_ranks(V2(:,2));

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
plot(ranks_from_dmaps,100*V2(:,3),'.')
hold on
%draw_tol = 0.004;
im_delta = 2;
hold on
for i=1:m
        image('cdata', image_set_aligned(:, :, :, i), 'xdata', [ranks_from_dmaps(i) ranks_from_dmaps(i)+im_delta], 'ydata', [100*V2(i,3) 100*V2(i,3)+im_delta])
end
axis equal

%% VDM

eps = median(W(:));
neigs = 12;

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

figure;
bar(diag(D))

%%

idx1 = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);
idx2 = find(embed_idx(1,:) == 7 & embed_idx(2,:) == 3);

figure;
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.')

%%

figure;
plot(embed_coord(wt,idx1),embed_coord(wt,idx2),'.')
hold on
plot(embed_coord(mut,idx1),embed_coord(mut,idx2),'.r')

%%

ranks_from_vdm = compute_ranks(embed_coord(:,idx1));

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
plot(ranks_from_vdm,1000*embed_coord(:,idx2),'.')
hold on
dist_in_dmaps = squareform(pdist(embed_coord(:, [idx1 idx2])));
%draw_tol = 0.004;
im_delta = 10;
hold on
for i=1:m
        image('cdata', image_set_aligned(:, :, :, i), 'xdata', [ranks_from_vdm(i) ranks_from_vdm(i)+im_delta], 'ydata', [1000*embed_coord(i,idx2) 1000*embed_coord(i,idx2)+im_delta])
end
axis equal

%%
[~, I ] =sort(embed_coord(:,idx1));

figure;
for i=1:m
    im1 = image_set_aligned(:,:,:,I(i));
    
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
end
    



