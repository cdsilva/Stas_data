clear all; 
close all;

% make green colormap
% cm_green = gray;
% cm_green(:,1) = 0;
% cm_green(:,3) = 0;

%%
% directories where images are stored
%ndirs = 3;
%image_dir = {'../membrane'; '../membrane2'; '../membrane3'};
%length_date = {'mar15'; 'oct16'; 'feb11'}; 

ndirs = 1;
image_dir = {'../membrane2'};
length_date = {'oct16'}; 

% pictures are in image channel 2
image_channel = 2;
% number of pixels (subsample images)
npixels = 100;

% load in membrane length data
mem_lengths = [];
nimages_vec = [];
for i=1:ndirs
    load(sprintf('%s/length_%s.mat', image_dir{i}, length_date{i}));
    mem_lengths = [mem_lengths; length(:,2)];
    nimages_vec = [nimages_vec; size(length, 1)];
    clear length;
end

% number of images
nimages = sum(nimages_vec);

%set image plotting parameters
subplot_dim1 = ceil(sqrt(nimages));
subplot_dim2 = ceil(nimages / subplot_dim1);

% store images in 3d array
image_set = zeros(npixels, npixels, nimages);

% load in images
figure;
for j=1:ndirs
    for i=1:nimages_vec(j)
        % read image
        im1 = imread(sprintf('%s/dpERK_staining/emb%02d.tif', image_dir{j}, i));
        
        % resize image
        im1 = imresize(im1, [npixels npixels]);
        
        % extract relevent color from image
        im1 = im1(:,:,image_channel);
        
        
        % index to store image
        idx = sum(nimages_vec(1:j-1))+i;
        
        % plot image
        subplot(subplot_dim1,subplot_dim2,idx)
        imshow(im1)
        
        %store image
        im1 = double(im1);
        image_set(:,:,idx) = im1;
    end
end

image_set(:,:,4) = image_set(:,:,3);
mem_lengths(4) = mem_lengths(3);

% save('image_data.mat', 'image_set', 'nimages', 'mem_lengths', 'npixels');
% 
% %%
% load image_data.mat

%%

% compute scattering transform of first image
x = image_set(:,:,1);
Wop = wavelet_factory_3d(size(x));
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% compute scattering transform of first image rotated
x2 = imrotate(x, 45, 'crop');
Wop2 = wavelet_factory_3d(size(x2));
Sx2 = scat(x2, Wop2);
Sx_mat2 = format_scat(Sx2);

%mean(mean(format_scat(Sx),2),3) - mean(mean(format_scat(Sx2),2),3);

%%

% set scattering transform parameters
filt_opt = struct();
filt_rot_opt = struct();
% oversampling is set to infty
% so that scattering of original and rotated
% image will be sampled at exactly the same points
scat_opt = struct();
scat_opt.oversampling = 10;

% store scattering invariants (one row per image)
sx_all = zeros(nimages, size(Sx_mat, 1));

% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);

% compute scattering invariants for each image
for i=1:nimages
    x = image_set(:,:,i);
    
    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%%

% dmaps

W = squareform(pdist(sx_all)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

figure;
plot(mem_lengths, V(:,2),'.')
xlabel('membrane length')
ylabel('\phi_2')

figure;
[~, I] = sort(V(:,2));
for i=1:nimages
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set(:,:,I(i))))
end

%%
W = squareform(pdist(sx_all, 'seuclidean')).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

figure;
plot(mem_lengths, V(:,2),'.')
xlabel('membrane length')
ylabel('\phi_2')

figure;
[~, I] = sort(V(:,2));
for i=1:nimages
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set(:,:,I(i))))
end



