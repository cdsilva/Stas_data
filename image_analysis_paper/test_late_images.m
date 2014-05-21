clear all
close all

%% directories where things are stored

im_save_dir = 'paper_figures2';

image_dir = '../membrane_pictures/14_0501_dpERK_late';

m = 132;
ind = 1:m;

%% image parameters

npixels = 100;

%set image plotting parameters
subplot_dim1 = floor(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim2:-1:1)/subplot_dim2, (1:subplot_dim1)/subplot_dim1);

nchannels = 3;
image_set = zeros(npixels, npixels, nchannels, m, 'uint8');
image_set_raw = zeros(1024, 1024, nchannels, m, 'uint8');

%% load images
figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    image_set_raw(:, :, :, i) = im1;
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    %     blur image
    H = fspecial('disk',5);
    im1(:, :, 1) = imfilter(im1(:, :, 1),H,'replicate');
    im1(:, :, 1) = 0.5*imadjust(im1(:, :, 1));
    
    im1(:, :, 3) = 1.5 * im1(:, :, 3);
    
    im1 = circshift(im1,[0 0 -1]);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    %subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    %store image
    image_set(:, :, :, i) = im1;
    
end

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 10;
shift_step = 2;
dim = 3;

% matlabpool open 4;
% [R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% matlabpool close;
%
%save('pairwise_alignments_late_noblur.mat', 'R', 'W');
%save('pairwise_alignments_late.mat', 'R', 'W');
%save('pairwise_alignments_nonuclei.mat', 'R', 'W');



load('pairwise_alignments_late.mat');

%load('pairwise_alignments_nonuclei.mat');


%% select which images to use
% all
ind = 1:m;

%ind = setdiff(1:m, [1 2 3 7 17 29]);
% used for pictures on 5/6
%ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 50 102]);

% use for pictures
%ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 32 102]);
%%THIS ONE
ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105 54 131 32 102 101 100]);

% used for pictures on 5/5
%ind = setdiff(1:m, [1 2 3 5 12 17 19 21 29 35 57 64 71 74 76 92 93 105]);

% used for pictures on 5/7
%ind = setdiff(1:m, [1 2 3 35 57 64 71 74 76 92 93 105]);

% only young embryos
%ind = [55,29,123,35,85,19,54,58,93,69,74,17,105,7,70,14,131,115,102,67,82,112,97,117,81,40,127,63,124,96,104,48,120,121,5,92,111,41,34,1,3,2,9,130,20,21,83,108,76,60,44,62,71,30,26,49,36,64,12,57];

image_set = image_set(:, :, :, ind);
image_set_raw = image_set_raw(:, :, :, ind);

W = W(ind, ind);

R_ind = reshape(dim*repmat(ind, dim, 1) - repmat((dim-1:-1:0)', 1, length(ind)), [], 1)';
R = R(R_ind, R_ind);

m = length(ind);

% subplot_dim1 = floor(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);
subplot_dim1 = 18;
subplot_dim2 = 6;

[Y, X] = meshgrid((subplot_dim2:-1:1)/subplot_dim2, (1:subplot_dim1)/subplot_dim1);

%% synchronization

R_opt = ang_synch(R, dim);

%%

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    im1 = image_set(:,:,:,i);
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    %subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
end
saveas(gcf,sprintf('%s/raw_data2', im_save_dir), 'pdf')


image_set_aligned = zeros(size(image_set), 'uint8');
image_set_raw_aligned = zeros(size(image_set_raw), 'uint8');
for i=1:m
    im_tmp = image_set(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set_raw(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_raw_aligned(:,:,:,i) = im1;
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    im1 = image_set_raw_aligned(:,:,:,i);
    
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    imshow(im1);
    
end

%% order

W2 = zeros(m);
for i=1:m
    imagei = image_set_aligned(:,:,:,i);
    for j=1:i-1
        imagej = image_set_aligned(:,:,:,j);
        d2 = sum((double(imagei(:))-double(imagej(:))).^2);
        W2(i, j) = d2;
        W2(j, i) = W2(i, j);
    end
end

eps2 = median(W2(:))/10;
neigs = 5;

[V, D] = dmaps(W, eps2, neigs);

[~, I] = sort(V(:,2));
figure;
set(gcf, 'paperposition', [0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1 1/subplot_dim2])
    
    imshow(image_set_raw_aligned(:,:,:,I(i)));
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1 1/subplot_dim2])
    
    imshow(image_set_aligned(:,:,:,I(i)));
end

%% try PCA

[eigenimages, D_PCA, proj_coeff] = PCA_images(image_set_aligned, m);

[~, I] = sort(proj_coeff(:, 1));
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    im1 = image_set_aligned(:,:,:,I(i));
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    %subplot(subplot_dim1, subplot_dim2, i);
    %subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(im1);
    
end
saveas(gcf,sprintf('%s/data2_PCA_ordered', im_save_dir), 'pdf')


figure;
for i=1:9
    subplot(3, 3, i)
    imshow(eigenimages(:, :, :, i))
end

nstages = 12;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
    make_subplot(nstages, 1, 0.01, i);
    stage_indices = I(max(1, round((i-1)*m/nstages)+1):min(m,round(i*m/nstages)));
    im1 = uint8(mean(double(image_set_aligned(:,:,:,stage_indices)), 4));
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 4])
set(gcf, 'paperposition',[0 0 8 4])
bar(diag(D_PCA(1:10, 1:10)) / sum(diag(D_PCA)))
xlabel('k')
ylabel('\lambda_k')
saveas(gcf,sprintf('%s/data2_PCA_variance', im_save_dir), 'pdf')

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 4])
set(gcf, 'paperposition',[0 0 8 4])
plot(proj_coeff(:,1),proj_coeff(:,2),'.')
xlabel('PCA projection 1')
ylabel('PCA projection 2')
saveas(gcf,sprintf('%s/data2_PCA_proj', im_save_dir), 'pdf')

%% VDM

eps = median(W(:))/10;
neigs = 6;
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
plot(diag(D),'.')
xlabel('k')
ylabel('\lambda_k')
saveas(gcf,sprintf('%s/data2_evals', im_save_dir), 'pdf')


figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
scatter(embed_idx(1,:), embed_idx(2, :), 500, var(embed_coord),'.')
colorbar
xlabel('k')
ylabel('l')
saveas(gcf,sprintf('%s/data2_coord_var', im_save_dir), 'pdf')


image_set_aligned = zeros(size(image_set), 'uint8');
image_set_raw_aligned = zeros(size(image_set_raw), 'uint8');
theta_adjust = 90;

for i=1:m
    im_tmp = image_set(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    im1 = circshift(im1, [5 10 0]);
    image_set_aligned(:,:,:,i) = im1;
    
    im_tmp = image_set_raw(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    im1 = imrotate(im1, theta_adjust, 'crop');
    im1 = circshift(im1, [50 100 0]);
    
    image_set_raw_aligned(:,:,:,i) = im1;
end

%% order using vdm

idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);

[~, I] = sort(embed_coord(:, idx));
if find(I == 10) < m/2
    I = flipud(I);
end
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    
    imshow(image_set_raw_aligned(:,:,:,I(i)));
end

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8*subplot_dim2/subplot_dim1])
set(gcf, 'paperposition',[0 0 8 8*subplot_dim2/subplot_dim1])
for i=1:m
    %subplot(subplot_dim1, subplot_dim2, i);
    %subplot('position', [X(i)-1/subplot_dim1 Y(i)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01])
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    im1 = image_set_aligned(:,:,:,I(i));
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    
    imshow(im1);
end
%saveas(gcf,sprintf('%s/VDM_ordered', im_save_dir), 'pdf')

%%

% nstages = 12;
%
% figure;
% set(gcf, 'paperunits', 'centimeters')
% set(gcf, 'papersize', [8 8/nstages])
% set(gcf, 'paperposition',[0 0 8 8/nstages])
% for i=1:nstages
%
%     %subplot('position', [(i-1)/nstages 0 1/nstages-0.005 1])
%     make_subplot(nstages, 1, 0.01, i);
%     stage_indices = I(max(1, round((i-1)*m/nstages)+1):min(m,round(i*m/nstages)));
%     im1 = uint8(mean(double(image_set_aligned(:,:,:,stage_indices)), 4));
%     im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
%     im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
%
%     imshow(im1,'initialmagnification','fit','border','tight')
% end
% saveas(gcf,sprintf('%s/average_trajectory', im_save_dir), 'pdf')


%%

frame_points = linspace(1, m, nstages);
window_eps = 20;
window_tol = 0.01;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8/nstages])
set(gcf, 'paperposition',[0 0 8 8/nstages])
for i=1:nstages
    
    window_weights = exp(-(frame_points(i)-(1:m)).^2/window_eps);
    window_weights = window_weights / sum(window_weights);
    window_weights(window_weights < window_tol) = 0;
    im1 = window_weights(1) * double(image_set_raw_aligned(:,:,:,I(1)));
    for j=2:m
        if window_weights(j) > window_tol
            im1 = im1 + window_weights(j) * double(image_set_raw_aligned(:,:,:,I(j)));
        end
    end
    im1 = uint8(im1 / sum(window_weights));
    
    im1 = circshift(im1,[0 0 -1]);
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    
    make_subplot(nstages, 1, 0.01, i);
    imshow(im1,'initialmagnification','fit','border','tight')
    
end
saveas(gcf,sprintf('%s/average_trajectory', im_save_dir), 'pdf')

%%

writerObj = VideoWriter('gastrulation.avi');
writerObj.FrameRate = 100;
open(writerObj);

nframes = 1000;
frame_points = linspace(1, m, nframes);
window_eps = 20;
window_tol = 0.01;

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8])

i = 1;
window_weights = exp(-(frame_points(i)-(1:m)).^2/window_eps);
window_weights = window_weights / sum(window_weights);
window_weights(window_weights < window_tol) = 0;
im1 = window_weights(1) * double(image_set_raw_aligned(:,:,:,I(1)));
for j=2:m
    if window_weights(j) > window_tol
        im1 = im1 + window_weights(j) * double(image_set_raw_aligned(:,:,:,I(j)));
    end
end
im1 = uint8(im1 / sum(window_weights));

im1 = circshift(im1,[0 0 -1]);
im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);

imshow(im1,'initialmagnification','fit','border','tight')
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

clf

for i=1:nframes
    
    window_weights = exp(-(frame_points(i)-(1:m)).^2/window_eps);
    window_weights = window_weights / sum(window_weights);
    window_weights(window_weights < window_tol) = 0;
    im1 = window_weights(1) * double(image_set_raw_aligned(:,:,:,I(1)));
    for j=2:m
        if window_weights(j) > window_tol
            im1 = im1 + window_weights(j) * double(image_set_raw_aligned(:,:,:,I(j)));
        end
    end
    im1 = uint8(im1 / sum(window_weights));
    
    im1 = circshift(im1,[0 0 -1]);
    im1(:, :, 1) = im1(:, :, 1)  + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2)  + im1(:, :, 3);
    
    imshow(im1,'initialmagnification','fit','border','tight')
    
    %pause(0.1)
    frame = getframe;
    writeVideo(writerObj,frame);
    clf
end

close(writerObj);
