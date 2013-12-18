% read in membrane pictures; align using vector dmaps

clear all
close all

res = '-r300';
fmt = '-djpeg';

cm_green = gray;
cm_green(:,1) = 0;
cm_green(:,3) = 0;

ndir = 3;
image_dir = {'../membrane','../membrane2','../membrane3'};
length_file = {'length_mar15','length_oct16','length_feb11'};
adjust_image = true;
image_channel = 2;
im_scale = 0.2;
max_deg = 20;

% dimension of rotations
rot_dim = 3;


%% load in membrane lengths
nimages = zeros(ndir, 1);
mem_lengths = [];
for j=1:ndir
    % load in membrane length data
    load(sprintf('%s/%s.mat', image_dir{j}, length_file{j}))
    mem_lengths = [mem_lengths; length(:,2)];
    nimages(j) = size(length, 1);
    clear length;
end

% where to store fourier expansion coefficients
ncoeff = sum(2*(0:max_deg)+1);
F_l = zeros(sum(nimages), ncoeff);
F_norm = zeros(sum(nimages), max_deg+1);
image_set = cell(sum(nimages),1);

%% load in images
for j=1:ndir
    for i=1:nimages(j)
        im1 = imread(sprintf('%s/emb%02d.tif', image_dir{j}, i));
        im1 = imresize(im1, im_scale);
        im1 = im1(:,:,image_channel);
        if adjust_image
            im1 = imadjust(im1);
            im1 = adjust_image_radially(im1);
        end
        im1 = double(im1);
        image_set{sum(nimages(1:j-1))+i} = im1;

        temp_F_l = zeros(1, ncoeff);

        for l=0:max_deg
            start_idx = sum(2*(0:(l-1))+1)+1;
            temp_F_l(start_idx:start_idx+2*l) = f_hat(l, im1);
            F_norm(sum(nimages(1:j-1))+i,l+1) = norm(temp_F_l(start_idx:start_idx+2*l));
        end

        F_l(sum(nimages(1:j-1))+i,:) = temp_F_l;
    end
end

%% dmaps

W = squareform(pdist(F_norm)).^2;

eps = median(median(W));

[V, D] = dmaps(W, eps, 10);

figure;
plot(V(:,2),mem_lengths,'.')

figure;
scatter(V(:,2),V(:,3),50,mem_lengths, '.')

delta1 = 0.04;
delta2 = 2;
figure;
plot(V(:,2),mem_lengths,'.')
hold on
for i=1:sum(nimages)
    imagesc([V(i,2) V(i,2)+delta1], [mem_lengths(i) mem_lengths(i)+delta2], image_set{i})
    colormap(gray)
end


