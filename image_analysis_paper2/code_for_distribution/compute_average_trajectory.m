function avg_images = compute_average_trajectory(images, nstages, window_scale)
% avg_images = compute_average_trajectory(images, nstages, window_scale)
% compute the average trajectory over a set of registered
% and ordered images
% images is a set of images, where images(:,:,:,i) is the i^th image
% nstages is the number of images in the average trajectory
% window_scale the kernel scale for computing the average images
% avg_images is the average trajectory, where avg_images(:,:,:,i) is the
% i^th image in the average trajectory

% check that images contains RGB images
if ndims(images) ~= 4
    disp('ERROR: average trajectory can only be computed over RGB images');
    return
end

[npixels, npixels2, nchannels, nimages] = size(images);

avg_images = zeros(npixels, npixels2, nchannels, nstages, 'like', images);

% points at which to compute the average images, relative to the true
% images
frame_points = linspace(1, nimages, nstages);

for i=1:nstages
    
    % compute weights for averaging images, using a Gaussian kernel
    window_weights = exp(-((1:nimages)-frame_points(i)).^2/window_scale^2);
    window_weights = window_weights / sum(window_weights);
    
    % add images with relevant weights
    im_tmp = immultiply(images(:,:,:,1), window_weights(1));
    for j=2:nimages
        im_tmp = imlincomb(1, im_tmp, window_weights(j), images(:,:,:,j));
    end
    
    % store average image
    avg_images(:,:,:,i) = im_tmp;
    
end