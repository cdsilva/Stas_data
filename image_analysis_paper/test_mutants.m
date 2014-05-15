clear all
close all

%% directories where things are stored

im_save_dir = 'paper_figures2';

% image_dir = '../membrane_pictures/snanull_big';
% m = 98;

% image_dir = '../membrane_pictures/snanull_old';
% m = 54;

image_dir = '../membrane_pictures/snanull_old2';
m = 47;
wt = [3 4 5 8  10 11 13 15 17 18 19 20 22 24 29 34 35 36 39 40 41 43 47]; %23
mut = [1 2 6 7 12 14 16 21 23 25 26 28 30 31 32 33 37 38 42 44 45 46]; %22

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

% matlabpool open 4;
% [R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% matlabpool close

%save('pairwise_alignments_snanull_old2.mat', 'R', 'W');
load('pairwise_alignments_snanull_old2.mat');


%% select which images to use

% ind = setdiff(1:m, [17]);
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

% W2 = W;
% for i=1:m
%     W2(i, i) = 0;
% end
%
% eps = median(W2(:));
%
% [V2, D] = dmaps(W2, eps, 10);
% %[V, D] = dmaps_weight(W,1,eps,10);
%
% figure;
% plot(V2(:,2),V2(:,3),'.')
%
% figure;
% im_delta = 0.05;
% plot(V2(:,2),V2(:,3),'.')
% hold on
% for i=1:m
%     image('cdata', image_set(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
%     hold on
% end
% axis equal
%
% return

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
%saveas(gcf,sprintf('%s/mutants_DMAPS_large', im_save_dir), 'pdf')

%% PCA after synch

[eigenimages, D_PCA, proj_coeff] = PCA_images(image_set_aligned, m);

diag(D_PCA) / sum(diag(D_PCA))


figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
im_delta = 300;
plot(proj_coeff(:,1),proj_coeff(:,2),'.')
dist_in_dmaps = squareform(pdist(proj_coeff(:, 1:2)));
draw_tol = 300;
hold on
for i=1:m
    if i == 1 || min(dist_in_dmaps(i, 1:i-1)) > draw_tol
        image('cdata', image_set_aligned(:, :, :, i), 'xdata', [proj_coeff(i,1) proj_coeff(i,1)+im_delta], 'ydata', [proj_coeff(i,2) proj_coeff(i,2)+im_delta])
        %hold on
   end
end
axis equal
saveas(gcf,sprintf('%s/mutants_PCA', im_save_dir), 'pdf')

%% VDM

eps = median(W(:));
neigs = 9;

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

%%

idx1 = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);
idx2 = find(embed_idx(1,:) == 7 & embed_idx(2,:) == 3);

if mean(embed_coord(:, idx1)) < 0
    embed_coord(:, idx1) = -embed_coord(:, idx1);
end
if mean(embed_coord(:, idx2)) > 0
    embed_coord(:, idx2) = -embed_coord(:, idx2);
end

figure;
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.')


[cluster_idx, C, SUMD] = kmeans(embed_coord(:, [idx1 idx2]), 3);
figure;
scatter(embed_coord(:,idx1),embed_coord(:,idx2),50, cluster_idx, '.')

init_idx = find(cluster_idx == find(C(:,1)<0 & C(:,2)<0));
mutant_idx = find(cluster_idx == find(C(:,1)>0));
wt_idx = find(cluster_idx == find(C(:,2)>0));

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.')
hold on
dist_in_dmaps = squareform(pdist(embed_coord(:, [idx1 idx2])));
draw_tol = 0.004;
im_delta = 0.004;
hold on
for i=1:m
    if i == 1 || min(dist_in_dmaps(i, 1:i-1)) > draw_tol
        image('cdata', image_set_aligned(:, :, :, i), 'xdata', [embed_coord(i,idx1) embed_coord(i,idx1)+im_delta], 'ydata', [embed_coord(i,idx2) embed_coord(i,idx2)+im_delta])
    end
end
axis equal
%saveas(gcf,sprintf('%s/mutants_VDM', im_save_dir), 'pdf')


% nstages = 5;
% figure;
% for j=1:3
%     curr_idx = find(cluster_idx == j);
%     m2 = length(curr_idx);
%     [~, I] = sort(embed_coord(curr_idx,idx2));
%     image_set_tmp = image_set_aligned(:, :, :, curr_idx);
%     set(gcf, 'paperunits', 'centimeters')
%     set(gcf, 'papersize', [8 8/nstages])
%     set(gcf, 'paperposition',[0 0 8 8/nstages])
%     for i=1:nstages
%         
%         %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
%         make_subplot(nstages, 3, 0.01, nstages*(j-1)+i);
%         stage_indices = I(max(1, round((i-1)*m2/nstages)+1):min(m2,round(i*m2/nstages)));
%         im1 = uint8(mean(double(image_set_tmp(:,:,:,stage_indices)), 4));
%         imshow(im1,'initialmagnification','fit','border','tight')
%     end
% end

nstages = 7;
figure;
curr_idx = [init_idx; mutant_idx];
m2 = length(curr_idx);
[~, I] = sort(embed_coord(curr_idx,idx1));
image_set_tmp = image_set_aligned(:, :, :, curr_idx);
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
    make_subplot(nstages, 1, 0.01, i);
    stage_indices = I(max(1, round((i-1)*m2/nstages)+1):min(m2,round(i*m2/nstages)));
    im1 = uint8(mean(double(image_set_tmp(:,:,:,stage_indices)), 4));
    imshow(im1,'initialmagnification','fit','border','tight')
end

figure;
curr_idx = [init_idx; wt_idx];
m2 = length(curr_idx);
[~, I] = sort(embed_coord(curr_idx,idx2));
image_set_tmp = image_set_aligned(:, :, :, curr_idx);
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
    make_subplot(nstages, 1, 0.01, i);
    stage_indices = I(max(1, round((i-1)*m2/nstages)+1):min(m2,round(i*m2/nstages)));
    im1 = uint8(mean(double(image_set_tmp(:,:,:,stage_indices)), 4));
    imshow(im1,'initialmagnification','fit','border','tight')
end

figure;
plot(embed_coord(wt,idx1),embed_coord(wt,idx2),'.')
hold on
plot(embed_coord(mut,idx1),embed_coord(mut,idx2),'.r')
% 45 is wrong

p_wt = polyfit(embed_coord(wt,idx1),embed_coord(wt,idx2),2);
p_mut = polyfit(embed_coord(mut,idx1),embed_coord(mut,idx2),2);

x_wt = linspace(min(embed_coord(wt,idx1)), max(embed_coord(wt,idx1)), 100);
x_mut = linspace(min(embed_coord(mut,idx1)), max(embed_coord(mut,idx1)), 100);

figure;
plot(embed_coord(wt,idx1),embed_coord(wt,idx2),'.')
hold on
plot(embed_coord(mut,idx1),embed_coord(mut,idx2),'.r')
plot(x_wt, polyval(p_wt,x_wt), '-b')
plot(x_mut, polyval(p_mut,x_mut), '-r')


%%

% embed_coord2 = embed_coord;
% for i=1:size(embed_coord, 2)
%     for j=1:2
%         embed_coord2(:, i) = embed_coord2(:, i) * sqrt(D(embed_idx(j,i), embed_idx(j,i)));
%     end
% end
% W2 = squareform(pdist(embed_coord2)).^2;
% eps2 = median(W2(:));
% 
% [V2, D2] = dmaps(W2, eps2, 10);
% 
% figure;
% plot(V2(:,2),V2(:,3),'.')
% 
% figure;
% im_delta = 0.05;
% set(gcf, 'paperunits', 'centimeters')
% set(gcf, 'papersize', [8 8])
% set(gcf, 'paperposition',[0 0 8 8]);
% plot(V2(:,2),V2(:,3),'.')
% hold on
% for i=1:m
%     image('cdata', image_set_aligned(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
%     hold on
% end
% axis equal




