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
    %im1(:, :, 2) = uint8(0.5*double(im1(:, :, 3)));
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

% matlabpool open 2;
% [R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% matlabpool close

%save('pairwise_alignments_snanull.mat', 'R', 'W');
load('pairwise_alignments_snanull.mat');


%% select which images to use

% ind = setdiff(1:m, [17]);
% ind = setdiff(1:m, [73 37 34 29 15 21 87 62 17 55 96 18]);
% ind = setdiff(1:m, [37 34 29 15 21 2 4 6]);
ind = setdiff(1:m, [37 34 29 15 21 2 4 6]);
mutation = mutation(ind);

image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);
nuclei = nuclei(:, :, ind);

W = W(ind, ind);

R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
R = R(R_ind, R_ind);

m = length(ind);

subplot_dim1 = floor(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

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
saveas(gcf,sprintf('%s/raw_data3b', im_save_dir), 'pdf')

%% synchronization

% R_opt = ang_synch(R, dim);
%
% image_set_aligned = zeros(size(image_set), 'uint8');
% image_set_aligned_withnuclei = zeros(size(image_set), 'uint8');
% for i=1:m
%     im_tmp = image_set(:,:,:,i);
%     %         for j=1:3
%     %             im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
%     %         end
%     im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
%     image_set_aligned(:,:,:,i) = im1;
%
%     im_tmp = image_set(:,:,:,i);
%     for j=1:3
%         im_tmp(:,:,j) = im_tmp(:,:,j) + nuclei(:,:,i);
%     end
%     im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
%     image_set_aligned_withnuclei(:,:,:,i) = im1;
% end
%
% %% dmaps after synch
%
% figure;
% for i=1:m
%     make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
%
%     imshow(image_set_aligned(:,:,:,i));
% end
%
% W2 = zeros(m);
% for i=1:m
%     im1 = image_set_aligned(:,:,:,i);
%     for j=1:i-1
%         im2 = image_set_aligned(:,:,:,j);
%         W2(i, j) = sum((double(im1(:)) - double(im2(:))).^2);
%         W2(j, i) = W2(i, j);
%     end
% end
%
% eps2 = median(W2(:))*10;
%
% [V2, D] = dmaps(W2, eps2, 10);
%
% figure;
% plot(V2(:,2),V2(:,3),'.');
%
% figure;
% plot(V2(~mutation,2),V2(~mutation,3),'.');
% hold on
% plot(V2(mutation,2),V2(mutation,3),'.r');
%
% figure;
% im_delta = 0.04;
% plot(V2(:,2),V2(:,3),'.')
% hold on
% for i=1:m
%     image('cdata', image_set_aligned(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
%     hold on
% end
% axis equal
%
% figure;
% set(gcf, 'paperunits', 'centimeters')
% set(gcf, 'papersize', [8 8])
% set(gcf, 'paperposition',[0 0 8 8]);
% im_delta = 0.04;
% plot(V2(:,2),V2(:,3),'.')
% dist_in_dmaps = squareform(pdist(V2(:, 2:3)));
% draw_tol = 0.04;
% hold on
% for i=1:m
%     if i == 1 || min(dist_in_dmaps(i, 1:i-1)) > draw_tol
%         image('cdata', image_set_aligned(:, :, :, i), 'xdata', [V2(i,2) V2(i,2)+im_delta], 'ydata', [V2(i,3) V2(i,3)+im_delta])
%         %hold on
%     end
% end
% axis equal
%
%
% %%
% ranks_from_dmaps = compute_ranks(V2(:,2));
%
% figure;
% set(gcf, 'paperunits', 'centimeters')
% set(gcf, 'papersize', [8 8])
% set(gcf, 'paperposition',[0 0 8 8]);
% plot(ranks_from_dmaps,100*V2(:,3),'.')
% hold on
% %draw_tol = 0.004;
% im_delta = 2;
% hold on
% for i=1:m
%         image('cdata', image_set_aligned(:, :, :, i), 'xdata', [ranks_from_dmaps(i) ranks_from_dmaps(i)+im_delta], 'ydata', [100*V2(i,3) 100*V2(i,3)+im_delta])
% end
% axis equal

%% VDM

eps = median(W(:))/5;
neigs = 42;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

figure;
plot(diag(D),'.')

theta_adjust = -90;

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
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
plot(abs(diag(D)),'.')
xlabel('k')
ylabel('|\lambda_k|')
saveas(gcf,sprintf('%s/data3_evalsb', im_save_dir), 'pdf')

%%

idx1 = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);
idx2 = find(embed_idx(1,:) == 7 & embed_idx(2,:) == 1);

if mean(embed_coord(:, idx1)) < 0
    embed_coord(:, idx1) = -embed_coord(:, idx1);
end
if mean(embed_coord(:, idx2)) > 0
    embed_coord(:, idx2) = -embed_coord(:, idx2);
end

figure;
plot(embed_coord(:,idx1),embed_coord(:,idx2),'.')


figure;
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'x')
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'o')

%%
[~, I ] =sort(embed_coord(:,idx1));

figure;
for i=1:m
    im1 = image_set_aligned(:,:,:,I(i));
    
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
end

%%

mut_draw = [27 26 88 79 65];
wt_draw = [3 20 70 58 28];
nstages = 5;

wt_ind = find(embed_coord(:, idx2) > 0);
mut_ind = find(embed_coord(:, idx2) < 0);

figure;
plot(embed_coord(wt_ind,idx1),embed_coord(wt_ind,idx2),'.')
hold on
plot(embed_coord(mut_ind,idx1),embed_coord(mut_ind,idx2),'o')

p_wt = polyfit(embed_coord(wt_ind,idx1),embed_coord(wt_ind,idx2),1);
p_mut = polyfit(embed_coord(mut_ind,idx1),embed_coord(mut_ind,idx2),1);

x_wt = linspace(-0.002, 0.02, 100);
x_mut = linspace(-0.002, 0.055, 100);

figure;
plot(embed_coord(~mutation,idx1),embed_coord(~mutation,idx2),'x')
hold on
plot(embed_coord(mutation,idx1),embed_coord(mutation,idx2),'o')
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
curve_delta_wt = 0.01;
curve_delta_mut = 0.01;
% patch([x_wt fliplr(x_wt)], [polyval(p_wt,x_wt)+curve_delta_wt polyval(p_wt,fliplr(x_wt))-curve_delta_wt], 'b', 'facealpha', 0.5, 'edgecolor','none')
% patch([x_mut fliplr(x_mut)], [polyval(p_mut,x_mut)+curve_delta_mut polyval(p_mut,fliplr(x_mut))-curve_delta_mut], 'r', 'facealpha', 0.5, 'edgecolor','none')
axis([-0.002 0.055 -0.06 0.02])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('first VDM coordinate')
ylabel('second VDM coordinate')
lgnd = legend('wild type','mutant', 'location','northeast');
set(lgnd,'fontsize',6);
saveas(gcf,sprintf('%s/mut_wt_vdm_embeddingb', im_save_dir), 'pdf')

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
% curve_delta_wt = std(embed_coord(wt_ind,idx2)-polyval(p_wt, embed_coord(wt_ind,idx2)));
% curve_delta_mut = std(embed_coord(mut_ind,idx2)-polyval(p_mut, embed_coord(mut_ind,idx2)));
curve_delta_wt = 0.011;
curve_delta_mut = 0.011;
patch([x_wt fliplr(x_wt)], [polyval(p_wt,x_wt)+curve_delta_wt polyval(p_wt,fliplr(x_wt))-curve_delta_wt], 'b', 'facealpha', 0.5, 'edgecolor','none')
hold on
patch([x_mut fliplr(x_mut)], [polyval(p_mut,x_mut)+curve_delta_mut polyval(p_mut,fliplr(x_mut))-curve_delta_mut], 'r', 'facealpha', 0.5, 'edgecolor','none')
axis([-0.002 0.055 -0.06 0.02])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
% xlabel('first VDM coordinate')
% ylabel('second VDM coordinate')
% lgnd = legend('wild type','mutant', 'location','northeast');
% set(lgnd,'fontsize',6);
saveas(gcf,sprintf('%s/mut_wt_vdm_embedding_backgroundb', im_save_dir), 'pdf')


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
saveas(gcf,sprintf('%s/wt_trajectoryb', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    make_subplot(nstages, 1,  0.01, i);
    im1 = image_set_aligned_withnuclei(:,:,:,mut_draw(i));
    im1(:, :, 2) = im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
end
saveas(gcf,sprintf('%s/mut_trajectoryb', im_save_dir), 'pdf')



