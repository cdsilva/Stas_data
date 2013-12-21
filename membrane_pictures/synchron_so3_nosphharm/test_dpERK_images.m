clear all
close all

%% parameters
dpERK_data = '../../membrane_lengths/oct16.mat';
dpERK_image_dir = '../membrane2/dpERK_staining';
dpERK_membrane_dir = '../membrane2';

% pictures are in image channel 2
image_channel = 2;
% number of pixels (subsample images)
npixels = 100;

load(dpERK_data, 'L');
m = size(L, 1);

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

% amont of border or buffer to add around the images
buffer_size = 50;

%size of "portion" of sphere on which to project
angle_proj = pi/2;

% total number of pixels
n = npixels+2*buffer_size;

% dimension of rotations
dim = 3;

% store images in 3d array
image_set = zeros(n, n, m);

%% load in images
figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', dpERK_image_dir, i));

    % resize image
    im1 = imresize(im1, [npixels npixels]);

    % extract relevent color from image
    im1 = im1(:,:,image_channel);

    %store image
    im1 = double(im1);
    image_set(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels,i) = im1;
end

%% compute pairwise alignments
tic
[R, W, angles] = align_data_nosph(image_set, angle_proj);
toc

%% vdm
%R_opt = ang_synch(R, 3);
eps = median(W(:));
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, 5);

%% plot
image_set_aligned = zeros(size(image_set));
figure;
for i=1:m
    image_tmp = image_set(:,:,i);
    
    R_all = R_opt(1:3, 1:3)'*R_opt(3*i-2:3*i,:);
    
    [alpha, beta, gamma] = calc_angles(R_all);    
    image_tmp = shift_image(image_tmp, alpha, beta, gamma, angle_proj);
    image_set_aligned(:,:,i) = image_tmp;
    
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_tmp(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels)),'InitialMagnification', 'fit')
end

%%
W2 = zeros(m);
for i=1:m
    for j=1:i-1
        W2(i,j) = sum(sum((image_set_aligned(:,:,i) - image_set_aligned(:,:,j)).^2));
        W2(j,i) = W2(i,j);
    end
end
eps2 = median(W2(:));
[V, D] = dmaps(W2, eps2, 5);

figure;
plot(V(:,2),L(:,1),'.')

%%

[~, embed_ind] = max(abs(corr(embed_coord, L(:,1))));

figure;
plot(embed_coord(:,embed_ind), L(:,1),'.')


% figure; 
% for i=1:size(embed_coord, 2)
%     plot(embed_coord(:,i), L(:,1),'.')
%     pause
%     clf
% end

[~, idx] = sort(embed_coord(:,embed_ind));

figure;
for i=1:m
    
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_set_aligned(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels,idx(i))),'InitialMagnification', 'fit')
end
