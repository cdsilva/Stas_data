clear all
close all

%% set parameters 

time_data = '../membrane_lengths/oct16.mat';
image_dir = '../membrane_pictures/membrane2/dpERK_staining';

%% load membrane lengths

load(time_data);
mem_lengths = L(:,1);
clear length;

m = length(mem_lengths);
ind = 1:m;

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% image parameters

npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

image_set = zeros(npixels, npixels, m);


%% load images

image_channel = [2, 3];
mat_file_names = {'dpERK_aligned_small.mat'; 'Dl_aligned_small.mat'};
mat_file_names_all = {'dpERK_aligned_small_all.mat'; 'Dl_aligned_small_all.mat'};

for j=1:length(image_channel)
    for i=1:m
        % read image
        im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
        
        % resize image
        im1 = imresize(im1, [npixels npixels]);
        
        % extract relevent color from image
        im1 = im1(:,:,image_channel(j));
        
        
        %store image
        im1 = double(im1);
        image_set(:, :, i) = im1;
        
    end
    
    %% raw dpERK images-- synchronization
    
    angle_proj = pi/8;
    shift_max = 20;
    shift_step = 4;
    dim = 3;
    
    [R, W] = compute_pairwise_alignments(image_set, angle_proj, shift_max, shift_step);
    
    %% VDM
    
    eps = median(W(:));
    neigs = 5;
    
    [R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
    
    %%
    
    figure;
    for i=1:m
        im1 = uint8(image_set(:,:,i));
        subplot(subplot_dim1, subplot_dim2, i);
        imshow(im1);
    end
    
    image_set_aligned = zeros(size(image_set));
    for i=1:m
        im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
        image_set_aligned(:,:,i) = double(im1);
    end
    
    %% order using vdm
    
    idx = find(embed_idx(1,:) == 4 & embed_idx(2,:) == 1);
    
    if corr(mem_lengths, embed_coord(:, idx)) < 0
        embed_coord(:, idx) = -embed_coord(:, idx);
    end
    
    [~, I] = sort(embed_coord(:, idx));
    figure;
    for i=1:m
        subplot(subplot_dim1, subplot_dim2, i);
        imshow(uint8(image_set_aligned(:,:,I(i))));
    end
    
    ranks_from_membranes = compute_ranks(mem_lengths);
    ranks_from_vdm = compute_ranks(embed_coord(:,idx));
    
    
    figure;
    plot(ranks_from_membranes, ranks_from_vdm, '.')
    xlabel('rank from membrane lengths')
    ylabel('rank from vdm')
    
    image_set_aligned2 = uint8(image_set_aligned(:,:,I));
    save(mat_file_names{j}, 'image_set_aligned2', 'ranks_from_membranes', 'ranks_from_vdm');
    save(mat_file_names_all{j})
end