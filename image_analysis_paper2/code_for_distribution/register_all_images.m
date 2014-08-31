function image_set2 = register_all_images(image_set, R_opt)
% image_set2 = register_all_images(image_set, R_opt)
% register all images using rotation matrices
% image_set are the images to register, where image_set(:,:,:,i)
% (for RGB images) or image_set(:,:,i) (for grayscale images) is the i^th
% image in the data set
% R_opt contains the optimal rotation matrices, returned by the vdm
% function
% image_set2 is the images from image_set, registered using the rotations in R_opt

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

