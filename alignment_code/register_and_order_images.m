function image_set_ordered = register_and_order_images(image_set)
% image_set is an npixels x npixels x nchannels x nimages array
% image_set_ordered is the same data set, registered and ordered using
% vector diffusion maps

[npixels, ~, nchannels, m] = size(image_set);

%% compute pairwise alignments

angle_proj = pi/8;
shift_max = 10;
shift_step = 2;
dim = 3;

% matlabpool open 2;
[R, W] = compute_pairwise_alignments_color(image_set, angle_proj, shift_max, shift_step);
% matlabpool close;

%% VDM

eps = median(W(:))/10;
neigs = 3*m;
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

image_set_aligned = zeros(size(image_set), 'uint8');

for i=1:m
    im_tmp = image_set(:,:,:,i);
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), im_tmp, angle_proj);
    image_set_aligned(:,:,:,i) = im1;
end

%% order using vdm

idx = 0;
var_tmp = 0;
for i=1:size(embed_idx, 2)
    if var(embed_coord(:,i)) > var_tmp && ...
            (embed_idx(1,i) == 4 || embed_idx(1,i) == 5 || embed_idx(1,i) == 6) && ...
            (embed_idx(2,i) == 1 || embed_idx(2,i) == 2 || embed_idx(2,i) == 3)
        idx = i;
        var_tmp = var(embed_coord(:,i));
    end
end

[~, I] = sort(embed_coord(:, idx));
image_set_ordered = image_set_aligned(:, :, :, I);

