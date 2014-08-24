function image_set2 = register_all_images(image_set, R_opt)

dim = size(R_opt, 2);

ndims_image_set = ndims(image_set);
nimages = size(image_set, ndims_image_set);

if size(R_opt, 1) ~= dim*nimages
    disp('ERROR: size of R_opt and image_set do not match');
    return
end

image_set2 = zeros(size(image_set), 'like', image_set);

if ndims_image_set == 4
    for i=1:nimages
        R_tmp = R_opt(dim*(i-1)+1:dim*i, :)';
        image_set2(:, :, :, i) = rotate_image(image_set(:,:, :,i), R_tmp);
    end
else % ndims_image_set == 3
    for i=1:nimages
        R_tmp = R_opt(dim*(i-1)+1:dim*i, :)';
        image_set2(:, :, i) = rotate_image(image_set(:,:,i), R_tmp);
    end
end

