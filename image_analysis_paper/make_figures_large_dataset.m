clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';
print_figures = false;

dpERK_data = '../membrane_pictures/large_dataset/time.mat';
dpERK_image_dir = '../membrane_pictures/large_dataset';
%dpERK_membrane_dir = '../membrane_pictures/membrane2';

%% load membrane lengths

load(dpERK_data);
mem_lengths = length;
clear length;

m = length(mem_lengths);
m = 16;

%% load images

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);

image_channel = 2;

% image indices to plot
% im_save_idx = [1,19,25,30,35,40,45,49];

% figure;
% for i=1:m
%     % read image
%     im1 = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir, i));
% 
%     % resize image
%     im1 = imresize(im1, [npixels npixels]);
% 
%     % extract relevent color from image
%     im1 = im1(:,:,image_channel);
% 
%     subplot(subplot_dim1, subplot_dim2, i)
%     imshow(im1);
%     
%     %store image
%     im1 = double(im1);
%     image_set(:, :, i) = im1;
%     
% end

%%
[~, ind] = max(mem_lengths);

im_template = imread(sprintf('%s/lat%02d.tif', dpERK_image_dir,ind));
%resize image
im_template = imresize(im_template, [npixels npixels]);
%extract relevent color from image
im_template = im_template(:,:,image_channel);

theta_vec = linspace(0, 360, 21);
theta_vec = theta_vec(1:end-1);

thetas = zeros(m,1);
xdisp = zeros(m,1);
ydisp = zeros(m,1);

figure;
for i=1:m
    thetas(i) = theta_vec(randi(20,1));
    xdisp(i) = randi([-10 10], 1);
    ydisp(i) = randi([-10 10], 1);
    
    im1 = imrotate(im_template, thetas(i), 'crop');
    im1 = circshift(im1, [xdisp(i) ydisp(i)]);
    
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(im1);
    
    %store image
    im1 = double(im1);
    image_set(:, :, i) = im1;
    
end


%%
% [~, ind] = sort(mem_lengths);
% ind = ind(end-10:end);
% 
% mem_lengths = mem_lengths(ind);
% image_set = image_set(:,:,ind);
% m = length(mem_lengths);
% 
% %set image plotting parameters
% subplot_dim1 = ceil(sqrt(m));
% subplot_dim2 = ceil(m / subplot_dim1);
% 
% figure;
% for i=1:m
%     subplot(subplot_dim1, subplot_dim2, i)
%     imshow(uint8(image_set(:,:,i)));
% end
% 
% L = mem_lengths;

%% raw dpERK images-- synchronization

angle_proj = pi/8;
shift_max = 20;
shift_step = 4;
dim = 3;
[R, W] = compute_pairwise_alignments(image_set, angle_proj, shift_max, shift_step);

%%

eps = median(W(:));
%neigs = 5;
%[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
%R_opt = ang_synch(R, 3);

[V, D] = eigs(R, dim);
if det(V(1:3,:)*D*V(1:3,:)') < 0
    V(:,end) = -V(:,end);
end

%R_opt =  V*sqrt(D);
% R_opt = zeros(size(V));
% for i=1:m
%     R_opt(dim*(i-1)+1:dim*i,:) = V(dim*(i-1)+1:dim*i,:) * inv(sqrtm(V(dim*(i-1)+1:dim*i,:)'*V(dim*(i-1)+1:dim*i,:)));
%     det(R_opt(dim*(i-1)+1:dim*i,:))
% end

%R2 = V*D*V';

R_opt = V*D*V(1:3,:)';
% R_opt = zeros(size(V));
% for i=1:m
%     [u, s, v] = svd(V(dim*(i-1)+1:dim*i,:));
%     R_opt(dim*(i-1)+1:dim*i,:) = u * v';
% end
% 
% for i=1:m
%     [u, s, v] = svd(V(dim*(i-1)+1:dim*i,:));
%     R_opt(dim*(i-1)+1:dim*i,:) = R_opt(dim*(i-1)+1:dim*i,:) * R_opt(1:3,:);
% end

for i=1:m
    R_opt(dim*(i-1)+1:dim*i,:) = R_opt(dim*(i-1)+1:dim*i,:) /(det(R_opt(dim*(i-1)+1:dim*i,:))^(1/3));
end

%%


figure;
j = m-1;
for i=1:m
    %im1 = rotate_image(R_opt(1:dim,1:dim)'*R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
    det(R_opt(dim*(i-1)+1:dim*i,:))
    %det(R_opt(1:dim,1:dim)'*R_opt(dim*(i-1)+1:dim*i,:))
    %im1 = rotate_image(R2(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j), uint8(image_set(:,:,i)), angle_proj);
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(im1);
end


