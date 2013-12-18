clear all; 
close all;

% make green colormap
cm_green = gray;
cm_green(:,1) = 0;
cm_green(:,3) = 0;

%%
% directories where images are stored
% ndirs = 3;
% image_dir = {'../membrane'; '../membrane2'; '../membrane3'};
% length_date = {'mar15'; 'oct16'; 'feb11'}; 

ndirs = 1;
image_dir = {'../membrane2'};
length_date = {'oct16'}; 
% image_dir = {'../membrane'};
% length_date = {'mar15'}; 
% image_dir = {'../membrane3'};
% length_date = {'feb11'}; 

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
subplot_dim1 = floor(sqrt(nimages));
subplot_dim2 = ceil(nimages / subplot_dim1);

% store images in 3d array
image_set = zeros(npixels, npixels, nimages);

% load in images
figure;
for j=1:ndirs
    for i=1:nimages_vec(j)
        % read image
        im1 = imread(sprintf('%s/emb%02d.tif', image_dir{j}, i));
        
        % resize image
        im1 = imresize(im1, [npixels npixels]);
        
        % extract relevent color from image
        im1 = im1(:,:,image_channel);
        
        % adjust intensity of image
        %im1 = adjust_image_radially(im1);
        %im1 = imadjust(im1);
        if i==1 && j==1
            im1 = imadjust(im1);
            im_ref = im1;
        else
            im1 = imhistmatch(im1, im_ref);
        end
        
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

% save('image_data.mat', 'image_set', 'nimages', 'mem_lengths', 'npixels');

%%
% load image_data.mat
% 
% %set image plotting parameters
% subplot_dim1 = ceil(sqrt(nimages));
% subplot_dim2 = ceil(nimages / subplot_dim1);

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
neigs = 10;
[V, D] = dmaps(W, eps, neigs);

figure;
scatter(V(:,2),V(:,3),50,mem_lengths,'.')

regres_idx = 1:3;
B = regress(mem_lengths, V(:,regres_idx));

figure;
plot(V(:,regres_idx)*B, mem_lengths, '.')

[~, I] = sort(V(:,regres_idx)*B);
figure;
for i=1:nimages
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set(:,:,I(i))))
end

B2 = [B(1); -B(3); B(2)];
[~, I] = sort(V(:,regres_idx)*B2);
figure;
for i=1:nimages
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set(:,:,I(i))))
end

% for i=2:neigs
%     figure;
%     plot(mem_lengths, V(:,i),'.')
% end
% 
% for j=2:neigs
%     figure;
%     [~, I] = sort(V(:,j));
%     for i=1:nimages
%         subplot(subplot_dim1,subplot_dim2,i)
%         imshow(uint8(image_set(:,:,I(i))))
%     end
% end

return

W = squareform(pdist(sx_all, 'seuclidean')).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

for i=2:5
    figure;
    plot(mem_lengths, V(:,i),'.')
end

figure;
[~, I] = sort(V(:,2));
for i=1:nimages
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set(:,:,I(i))))
end



%%

% do PCA of scattering transform
nPCA = 10;

[V, D] = PCA(sx_all, nPCA);

a = (sx_all - repmat(mean(sx_all), nimages, 1)) * V;

figure;
plot(a(:,1), a(:,2), '.')

figure;
scatter(a(:,1), a(:,2), 50, mem_lengths, '.')

figure;
scatter3(a(:,1), a(:,2), a(:,3), 50, mem_lengths, '.')

% use regression to fit mem length function of PCA coeff
B = regress(mem_lengths, [ones(nimages, 1) a]);

figure;
plot(mem_lengths, [ones(nimages, 1) a] * B, '.')
xlabel('true membrange length')
ylabel('predicted membrane length')


%%

% do PCA of raw images

% reshape image array into 2d array (images are unraveled)
image_all = zeros(nimages, npixels^2);
for i=1:nimages
    image_all(i,:) = reshape(image_set(:,:,i), 1, []);
end

image_mean = mean(image_all);

[V, D] = PCA(image_all, 10);

a = (image_all - repmat(image_mean, nimages, 1)) * V;

B = regress(mem_lengths, [ones(nimages, 1) a]);

figure;
plot(mem_lengths, [ones(nimages, 1) a] * B, '.')
xlabel('true membrange length')
ylabel('predicted membrane length')
