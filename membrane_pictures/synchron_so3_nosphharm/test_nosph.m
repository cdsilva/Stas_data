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

% number of samples
m = 30;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

% maximum shift allowed in images
shift_max = 5;

buffer_size = 20;

%size of "portion" of sphere on which to project
angle_proj = pi/4;

% total number of pixels
n = npixels+2*buffer_size;

% dimension of rotations
dim = 3;

% store images in 3d array
image_set = zeros(n, n, m);

% store exact rotations
R_exact = zeros(dim*m, dim);

f1 = figure;
f2 = figure;
for i=1:m
    
    im_tmp = zeros(n);
    im_tmp(buffer_size+1:buffer_size+npixels, buffer_size+1:buffer_size+npixels) = double(im1);
    
    alpha = 360 * rand;
    beta = randi([-shift_max, shift_max]) * angle_proj / n;
    gamma = randi([-shift_max, shift_max]) * angle_proj / n;
    
    im_tmp = shift_image(im_tmp, alpha, beta, gamma, angle_proj);
    
    image_set(:,:,i) = im_tmp;
    
    figure(f1)
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(im_tmp))
    
    R_exact(3*i-2:3*i,:) = calc_rot_matrix(alpha, beta, gamma);
    
    R_all = R_exact(3*i-2:3*i,:)';
    [alpha, beta, gamma] = calc_angles(R_all);    
    
    figure(f2)
    
    im_tmp = shift_image(im_tmp, alpha, beta, gamma, angle_proj);
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(im_tmp))
    
end

%return

tic
[R, W, angles] = align_data_nosph(image_set, angle_proj);
toc

R_opt = ang_synch(R, 3);
%eps = median(W(:));
%[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, 5);

%%
image_set_aligned = zeros(size(image_set));

figure;
for i=1:m
    image_tmp = image_set(:,:,i);
    
    R_all = R_opt(1:3, 1:3)'*R_opt(3*i-2:3*i,:);
    %R_all = R_opt(3*i-2:3*i,:);
    %R_all = R_exact(3*i-2:3*i,:)';
    
    [alpha, beta, gamma] = calc_angles(R_all);    
    
    
    image_tmp = shift_image(image_tmp, alpha, beta, gamma, angle_proj);
    
    image_set_aligned(:,:,i) = image_tmp;
    
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_tmp))
end

figure;
for i=1:m
    imshow(uint8(image_set_aligned(:,:,i)))
    pause
end

return

%%
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
        title(sprintf('j=%d, norm = %2.2f', idx2, norm(image_tmp-image_set(:,:,idx2))))
        pause
    end
end



