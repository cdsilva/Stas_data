clear all
close all


dpERK_data = '../../membrane_lengths/oct16.mat';
dpERK_image_dir = '../membrane2/dpERK_staining';
dpERK_membrane_dir = '../membrane2';

% pictures are in image channel 2
image_channel = 2;
% number of pixels (subsample images)
npixels = 100;

im1 = imread(sprintf('%s/emb%02d.tif', dpERK_image_dir, 1));
% resize image
im1 = imresize(im1, [npixels npixels]);

% extract relevent color from image
im1 = im1(:,:,image_channel);

m = 10;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

shift_max = 4;
angle_proj = pi/2;

n = npixels+2*shift_max;
dim = 3;

% store images in 3d array
image_set = zeros(n, n, m);

R_exact = zeros(dim*m, dim);

f1 = figure;
f2 = figure;
for i=1:m
    
    im_tmp = zeros(n);
    im_tmp(shift_max+1:shift_max+npixels, shift_max+1:shift_max+npixels) = double(im1);
    
    rand_angle = 360*rand;
    
    im_tmp = imrotate(im_tmp, rand_angle, 'crop');
    
    trans = randi([-shift_max, shift_max], 1, 2);
    im_tmp = circshift(im_tmp, trans);
    
    image_set(:,:,i) = im_tmp;
    
    figure(f1)
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(im_tmp))
    
    alpha = rand_angle * pi / 180;
    beta = trans(1) / n * angle_proj;
    gamma = trans(2) / n * angle_proj;
    Rx = [1 0 0;
        0 cos(alpha) -sin(alpha);
        0 sin(alpha) cos(alpha)];
    
    
    Ry = [cos(beta) 0 sin(beta);
        0 1 0;
        -sin(beta) 0 cos(beta)];
    
    Rz = [cos(gamma) -sin(gamma) 0;
        sin(gamma) cos(gamma) 0;
        0  0 1];
    
    
    R_exact(3*i-2:3*i,:) = Rz * Ry * Rx;
    
    R_all = R_exact(3*i-2:3*i,:)';
        alpha = atan2(R_all(3,2), R_all(3,3));
    gamma = atan2(R_all(2,1), R_all(1,1));
    beta = -asin(R_all(3,1));
    
    figure(f2)
    %im_tmp = circshift(im_tmp, -trans);
    %im_tmp = imrotate(im_tmp, -rand_angle, 'crop');
    
    im_tmp = imrotate(im_tmp, alpha*(180/pi), 'crop');
    im_tmp = circshift(im_tmp, round([beta gamma] * (npixels+2*shift_max)/angle_proj));
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(im_tmp))
    
end

%return

% %store image
% im1 = double(im1);
% image_set(shift_max+1:shift_max+npixels,shift_max+1:shift_max+npixels,i) = im1;
%
%
%
%
% m = 52;
%
% load(dpERK_data, 'L');
% [~, idx] = sort(L(:,1), 'descend');
%
% % pictures are in image channel 2
% image_channel = 2;
% % number of pixels (subsample images)
% npixels = 100;
%
% %set image plotting parameters
% subplot_dim1 = ceil(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);
%
% shift_max = 8;
%
% % store images in 3d array
% image_set = zeros(npixels+2*shift_max, npixels+2*shift_max, m);
%
% % load in images
% figure;
% for i=1:m
%     % read image
%     im1 = imread(sprintf('%s/emb%02d.tif', dpERK_image_dir, idx(i)));
%
%     % resize image
%     im1 = imresize(im1, [npixels npixels]);
%
%     % extract relevent color from image
%     im1 = im1(:,:,image_channel);
%
%     %store image
%     im1 = double(im1);
%     image_set(shift_max+1:shift_max+npixels,shift_max+1:shift_max+npixels,i) = im1;
% end

tic
[R, W, angles] = align_data_nosph(image_set, angle_proj);
toc

R_opt = ang_synch(R, 3);
%eps = median(W(:));
%[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, 5);

figure;
for i=1:m
    image_tmp = image_set(:,:,i);
    
    %R_all = (R_opt(1:3, 1:3)'*R_opt(3*i-2:3*i,:))';
    %R_all = R_opt(3*i-2:3*i,:)';
    R_all = R_exact(3*i-2:3*i,:)';
    
    alpha = atan2(R_all(3,2), R_all(3,3));
    gamma = atan2(R_all(2,1), R_all(1,1));
    beta = -asin(R_all(3,1));
    
    image_tmp = imrotate(image_tmp, (180/pi)*alpha, 'crop');
    image_tmp = circshift(image_tmp, round([beta gamma]*(npixels+2*shift_max)/angle_proj));
    
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_tmp))
end

return


figure;
for idx1=1:m
    for idx2=1:m
        R_all = R(3*(idx1-1)+1:3*idx1, 3*(idx2-1)+1:3*idx2);
        alpha = atan2(R_all(3,2), R_all(3,3));
        gamma = atan2(R_all(2,1), R_all(1,1));
        beta = -asin(R_all(3,1));
        
        image_tmp = image_set(:,:,idx1);
        image_tmp = imrotate(image_tmp, (180/pi)*alpha, 'crop');
        image_tmp = circshift(image_tmp, round([beta gamma]*(npixels+2*shift_max)/angle_proj));
        subplot(2,2, 1);
        imshow(uint8(image_tmp))
        title(sprintf('i=%d', idx1))
        subplot(2,2, 2);
        imshow(uint8(image_set(:,:,idx2)))
        title(sprintf('j=%d', idx2))
        pause
    end
end



