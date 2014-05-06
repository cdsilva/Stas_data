clear all
close all

%% directories where things are stored
time_data = '../membrane_lengths/oct16.mat';
image_dir = '../membrane_pictures/membrane2/dpERK_staining';

%% load membrane lengths

load(time_data);
mem_lengths = L(:,1);

% remove some embryos
ind = setdiff(1:52, [9, 15, 22]);

m = length(ind);

mem_lengths = mem_lengths(ind);

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% image parameters

npixels = 101;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim1:-1:1)/subplot_dim1, (1:subplot_dim2)/subplot_dim2);

data = zeros(npixels, npixels, m);

%% load images

for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    im1 = im1(:, :, 2);
    
    %store image
    data(:, :, i) = im1;
    
end

%%

r_max = floor(npixels/2);
n_nbor = 5;
isrann = 0;

[ class, rot, corr, FBsPCA_data, timing ] = Initial_classification_noReflection(data, r_max, n_nbor, isrann );
%[ class, class_refl, rot, corr, FBsPCA_data, timing ] = Initial_classification(data, r_max, n_nbor, isrann );

%%

k = n_nbor;
flag = 0;
[ class_VDM, angle, VV_LP ] = VDM_noReflection(class, corr, rot, k, flag, n_nbor) ;

%%

coord = zeros(m, 1);
theta = zeros(m, 1);

for i=1:m
    coord(i) = VV_LP(i, 2) * conj(VV_LP(i, 1));
    theta(i) = atan2(imag(VV_LP(i, 1)), real(VV_LP(i, 1))) * (180/pi);
end

[~, I] = sort(real(coord));

figure;
for i=1:m
    subplot(subplot_dim1, subplot_dim2, i)
    %imshow(uint8(data(:, :, I(i))));
    imshow(imrotate(uint8(data(:, :, I(i))), theta(I(i)), 'crop'))
end
    