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


mutation = true(m, 1);
mutation(wt) = false;

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
    %im1(:, :, 2) = 0;
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

% matlabpool open 2;
% [R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% matlabpool close

%save('pairwise_alignments_snanull_old2.mat', 'R', 'W');
load('pairwise_alignments_snanull_old2.mat');


%% select which images to use

% ind = setdiff(1:m, [4 15 19 21 47]);
ind = setdiff(1:m, [16 27 34 4 14]);

image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);
nuclei = nuclei(:, :, ind);

mutation = mutation(ind);

W = W(ind, ind);

R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
R = R(R_ind, R_ind);

m = length(ind);

subplot_dim1 = floor(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

%% save data set
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [subplot_dim1 subplot_dim2])
set(gcf, 'paperposition',[0 0 subplot_dim1 subplot_dim2])
for i=1:m
    im1 = image_set(:,:,:,i);
    
    for j=1:3
        im1(:,:,j) = im1(:,:,j) + nuclei(:,:,i);
    end

    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
    
end
%saveas(gcf,sprintf('%s/raw_data3', im_save_dir), 'pdf')

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
plot(V2(~mutation,2),V2(~mutation,3),'.')
hold on
plot(V2(mutation,2),V2(mutation,3),'.r')

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

% [eigenimages, D_PCA, proj_coeff] = PCA_images(image_set_aligned, m);
% 
% diag(D_PCA) / sum(diag(D_PCA))
% 
% 
% figure;
% set(gcf, 'paperunits', 'centimeters')
% set(gcf, 'papersize', [8 8])
% set(gcf, 'paperposition',[0 0 8 8]);
% im_delta = 300;
% plot(proj_coeff(:,1),proj_coeff(:,2),'.')
% dist_in_dmaps = squareform(pdist(proj_coeff(:, 1:2)));
% draw_tol = 300;
% hold on
% for i=1:m
%     if i == 1 || min(dist_in_dmaps(i, 1:i-1)) > draw_tol
%         image('cdata', image_set_aligned(:, :, :, i), 'xdata', [proj_coeff(i,1) proj_coeff(i,1)+im_delta], 'ydata', [proj_coeff(i,2) proj_coeff(i,2)+im_delta])
%         %hold on
%    end
% end
% axis equal
%saveas(gcf,sprintf('%s/mutants_PCA', im_save_dir), 'pdf')

%% VDM

eps = median(W(:))/2;
neigs = 3*m;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

theta_adjust = -55;

% mut_draw = [ 37 32 1 22 27];
% wt_draw = [19 7 10 9 34];
% nstages = 5;
% mut_draw = [ 37  1 22 27];
% wt_draw = [19  10 9 34];
% nstages = 4;
mut_draw = [29 5 1 27];
wt_draw = [21  10 9 34];
nstages = 4;

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
figure; 
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
plot(abs(diag(D)),'.')
xlabel('k')
ylabel('|\lambda_k|')
% saveas(gcf,sprintf('%s/data3_evals', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
nbins = 25;
hist(diag(D), nbins)
set(gca, 'xlim', [-1 1])
set(gca, 'ylim', [0 45])
grid on
set(gca, 'ytick', 0:3:45)
set(gca, 'xtick', -1:0.2:1)
xlabel('\lambda')
ylabel('Pr(\lambda)')
% hold on
% [nelements, centers] = hist(diag(D(7:end, 7:end)), nbins);
% [alpha, rad] = fit_semicircle(centers, nelements);
% plot(centers, wigner_semicircle(centers, rad, alpha))
saveas(gcf,sprintf('%s/data3_evals_dist', im_save_dir), 'pdf')

return

figure;  
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
scatter(embed_idx(1,:), embed_idx(2, :), 500, var(embed_coord),'.')
colorbar
xlabel('k')
ylabel('l')
% saveas(gcf,sprintf('%s/data3_coord_var', im_save_dir), 'pdf')


%% find relevant VDM coordinates

idx1 = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);
idx2 = find(embed_idx(1,:) == 7 & embed_idx(2,:) == 3);
% idx2 = find(embed_idx(1,:) == 7 & embed_idx(2,:) == 1);

if mean(embed_coord(:, idx1)) < 0
    embed_coord(:, idx1) = -embed_coord(:, idx1);
end
if mean(embed_coord(:, idx2)) < 0
    embed_coord(:, idx2) = -embed_coord(:, idx2);
end

figure;
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.')


figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.k', 'markersize', 10)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
saveas(gcf,sprintf('%s/data3_embed_bw', im_save_dir), 'pdf')


figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'xb', 'markersize', 5)
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'or', 'markersize', 5)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
lgnd = legend('wild type','mutant', 'location','northeast');
set(lgnd,'fontsize',6);
% saveas(gcf,sprintf('%s/data3_embed', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'.b', 'markersize', 10)
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'.',  'color', [0 0.8 0], 'markersize', 10)
plot(embed_coord(wt_draw,idx1),embed_coord(wt_draw,idx2),'ob', 'markersize', 8)
hold on
plot(embed_coord(mut_draw,idx1),embed_coord(mut_draw,idx2),'o', 'color', [0 0.8 0], 'markersize', 8)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
lgnd = legend('wild type','mutant', 'location','northeast');
set(lgnd,'fontsize',6);
saveas(gcf,sprintf('%s/data3_embed_color', im_save_dir), 'pdf')


return

%% draw curve with select images



figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
coord2_scale = 0.6;
plot(embed_coord(~mutation,idx1),coord2_scale*embed_coord(~mutation,idx2),'xk', 'markersize', 5)
hold on
plot(embed_coord(mutation,idx1),coord2_scale*embed_coord(mutation,idx2),'ok', 'markersize', 5)
hold on
im_delta = 0.004;
draw_delta = 0.008;
hold on
plot(embed_coord(mut_draw,idx1),coord2_scale*embed_coord(mut_draw,idx2),'ob', 'markersize', 5, 'linewidth', 1)
plot(embed_coord(wt_draw,idx1),coord2_scale*embed_coord(wt_draw,idx2),'xb', 'markersize', 5, 'linewidth', 1)
for i=mut_draw
   image('cdata', image_set_aligned(:, :, :, i), 'xdata', [embed_coord(i,idx1)-sign(embed_coord(i,idx1)-0.01)*draw_delta-im_delta embed_coord(i,idx1)-sign(embed_coord(i,idx1)-0.01)*draw_delta+im_delta], 'ydata', [coord2_scale*embed_coord(i,idx2)-draw_delta-im_delta coord2_scale*embed_coord(i,idx2)-draw_delta+im_delta])
end
for i=wt_draw
   image('cdata', image_set_aligned(:, :, :, i), 'xdata', [embed_coord(i,idx1)-draw_delta-im_delta embed_coord(i,idx1)-draw_delta+im_delta], 'ydata', [coord2_scale*embed_coord(i,idx2)-im_delta coord2_scale*embed_coord(i,idx2)+im_delta])
end
axis equal
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
lgnd = legend('wild type','mutant', 'location','northeast');
set(lgnd,'fontsize',6);
%print(sprintf('%s/mut_wt_vdm_embedding2', im_save_dir), '-dpdf', '-r600')

%% clustering

nclusters = 3;

cluster_idx = kmeans(embed_coord(:, [idx1 idx2]), nclusters, 'start', [-0.01 -0.01; 0.01 0.04; 0.04 -0.01]);
figure;
scatter(embed_coord(:,idx1),embed_coord(:,idx2),50, cluster_idx, '.')

wt_ind = [];
mut_ind = [];

figure;
symb = '.ox';
%scatter(embed_coord(:,idx1),embed_coord(:,idx2),50, cluster_idx, '.')
hold on
for i=1:nclusters
%for i=2:nclusters
    tmp_idx = find(cluster_idx == i);
    %tmp_idx = find(cluster_idx == i | cluster_idx == 1);
    data = embed_coord(tmp_idx, [idx1 idx2]);
    mean_data = mean(data);
    data = data - repmat(mean_data, size(data, 1), 1);
    
    [V,D]=PCA(data,1);
    proj_data = (data * V(:, 1) * V(:, 1)' + repmat(mean_data, size(data, 1), 1));
    r = compute_ranks(data * V(:, 1));
    [~, I] = sort(data * V(:, 1));
    scatter(embed_coord(tmp_idx,idx1),embed_coord(tmp_idx,idx2),50, r/max(r), symb(i));
    %plot(proj_data(:,1), proj_data(:,2), '.');
    
    if i==1
        wt_ind = [wt_ind; tmp_idx(I)];
        mut_ind = [mut_ind; tmp_idx(I)];
    elseif i==2 
        wt_ind = [wt_ind; tmp_idx(I)];
    elseif i==3
        mut_ind = [mut_ind; tmp_idx(I)];
    end
end

% wt_ind = setdiff(wt_ind, 1);
% mut_ind = [mut_ind; 1];

figure;
plot(embed_coord(wt_ind,idx1),embed_coord(wt_ind,idx2),'.')
hold on
plot(embed_coord(mut_ind,idx1),embed_coord(mut_ind,idx2),'o')

p_wt = polyfit(embed_coord(wt_ind,idx1),embed_coord(wt_ind,idx2),1);
p_mut = polyfit(embed_coord(mut_ind,idx1),embed_coord(mut_ind,idx2),2);

x_wt = linspace(-0.025, 0.05, 100);
x_mut = linspace(-0.025, 0.055, 100);

figure;
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'.')
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'.r')
plot(x_wt, polyval(p_wt,x_wt), '-b')
plot(x_mut, polyval(p_mut,x_mut), '-r')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'xk', 'markersize', 5)
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'ok', 'markersize', 5)
plot(embed_coord(wt_draw,idx1),embed_coord(wt_draw,idx2),'xk', 'markersize', 5, 'linewidth', 2)
hold on
plot(embed_coord(mut_draw,idx1),embed_coord(mut_draw,idx2),'ok', 'markersize', 5, 'linewidth', 2)
% curve_delta_wt = std(embed_coord(wt_ind,idx2)-polyval(p_wt, embed_coord(wt_ind,idx2)));
% curve_delta_mut = std(embed_coord(mut_ind,idx2)-polyval(p_mut, embed_coord(mut_ind,idx2)));
curve_delta_wt = 0.05;
curve_delta_mut = 0.01;
patch([x_wt fliplr(x_wt)], [polyval(p_wt,x_wt)+curve_delta_wt polyval(p_wt,fliplr(x_wt))-curve_delta_wt], 'b', 'facealpha', 0.5, 'edgecolor','none')
patch([x_mut fliplr(x_mut)], [polyval(p_mut,x_mut)+curve_delta_mut polyval(p_mut,fliplr(x_mut))-curve_delta_mut], 'r', 'facealpha', 0.5, 'edgecolor','none')
axis([-0.025 0.055 -0.07 0.07])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
lgnd = legend('wild type','mutant', 'location','northeast');
set(lgnd,'fontsize',6);

%% draw all images

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

%% draw with curve
    
wt_idx = find(embed_coord(:, idx2) > 2.75*embed_coord(:,idx1)+0.005);
mut_6dx = find(embed_coord(:, idx2) < 0.01);
% p_wt = [    3.0684    0.04];
% p_mut = [  -30.6192    0.8206    0.0039];
p_wt = [3.8522    0.0288];
p_mut = [-45    1.6066    0.005];

x_wt = linspace(-0.025, 0.05, 100);
x_mut = linspace(-0.025, 0.055, 100);

figure;
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'.')
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'.r')
plot(x_wt, polyval(p_wt,x_wt), '-b')
plot(x_mut, polyval(p_mut,x_mut), '-r')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'xk', 'markersize', 5)
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'ok', 'markersize', 5)
plot(embed_coord(wt_draw,idx1),embed_coord(wt_draw,idx2),'xk', 'markersize', 5, 'linewidth', 2)
hold on
plot(embed_coord(mut_draw,idx1),embed_coord(mut_draw,idx2),'ok', 'markersize', 5, 'linewidth', 2)
%curve_delta_wt = 0.028;
%curve_delta_mut = 0.021;
%patch([x_wt fliplr(x_wt)], [polyval(p_wt,x_wt)+curve_delta_wt polyval(p_wt,fliplr(x_wt))-curve_delta_wt], 'b', 'facealpha', 0.5, 'edgecolor','none')
%patch([x_mut fliplr(x_mut)], [polyval(p_mut,x_mut)+curve_delta_mut polyval(p_mut,fliplr(x_mut))-curve_delta_mut], 'r', 'facealpha', 0.5, 'edgecolor','none')
% axis([-0.025 0.055 -0.07 0.07])
axis([-0.02 0.055 -0.04 0.07])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
lgnd = legend('wild type','mutant', 'location','northeast');
set(lgnd,'fontsize',6);
% saveas(gcf,sprintf('%s/mut_wt_vdm_embedding', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 6])
set(gcf, 'paperposition',[0 0 8 6])
% plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'xk', 'markersize', 5)
% hold on
% plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'ok', 'markersize', 5)
% plot(embed_coord(wt_draw,idx1),embed_coord(wt_draw,idx2),'xk', 'markersize', 5, 'linewidth', 2)
% hold on
% plot(embed_coord(mut_draw,idx1),embed_coord(mut_draw,idx2),'ok', 'markersize', 5, 'linewidth', 2)
curve_delta_wt = 0.028;
curve_delta_mut = 0.018;
patch([x_wt fliplr(x_wt)], [polyval(p_wt,x_wt)+curve_delta_wt polyval(p_wt,fliplr(x_wt))-curve_delta_wt], 'b', 'facealpha', 0.5, 'edgecolor','none')
hold on
patch([x_mut fliplr(x_mut)], [polyval(p_mut,x_mut)+curve_delta_mut polyval(p_mut,fliplr(x_mut))-curve_delta_mut], 'r', 'facealpha', 0.5, 'edgecolor','none')
% axis([-0.025 0.055 -0.07 0.07])
axis([-0.02 0.055 -0.04 0.07])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
%xlabel('first VDM coordinate')
%ylabel('second VDM coordinate')
%lgnd = legend('wild type','mutant', 'location','northeast');
%set(lgnd,'fontsize',6);
%saveas(gcf,sprintf('%s/mut_wt_vdm_embedding_background', im_save_dir), 'pdf')
% print(sprintf('%s/mut_wt_vdm_embedding_background', im_save_dir), '-dpdf', '-r600')


%% draw images for select points
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    make_subplot(nstages, 1,  0.01, i);
    im1 = image_set_aligned_withnuclei(:,:,:,wt_draw(i));
    im1(:, :, 2) = im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
saveas(gcf,sprintf('%s/wt_trajectory', im_save_dir), 'pdf')


figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    make_subplot( nstages,1,  0.01, i);
    im1 = image_set_aligned_withnuclei(:,:,:,mut_draw(i));
    im1(:, :, 2) = im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
saveas(gcf,sprintf('%s/mut_trajectory', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8/nstages 8])
set(gcf, 'paperposition',[0 0 8/nstages 8])
for i=1:nstages
    
    make_subplot(1, nstages,  0.01, i);
    im1 = image_set_aligned_withnuclei(:,:,:,wt_draw(i));
    im1(:, :, 2) = im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
% saveas(gcf,sprintf('%s/wt_trajectory2', im_save_dir), 'pdf')


figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8/nstages 8])
set(gcf, 'paperposition',[0 0 8/nstages 8])
for i=1:nstages
    
    make_subplot(1, nstages,  0.01, i);
    im1 = image_set_aligned_withnuclei(:,:,:,mut_draw(i));
    im1(:, :, 2) = im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
% saveas(gcf,sprintf('%s/mut_trajectory2', im_save_dir), 'pdf')