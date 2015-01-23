function images = apply_image_functions(images_raw, dim, channel_weight, channel_blur, channel_normalize, channel_mean_center, resize_image)

images = images_raw;

nimages = size(images_raw, ndims(images_raw));
if dim == 3
    nstack = size(images_raw, ndims(images_raw)-1);
    if ndims(images_raw) == 5
        nchannels = size(images_raw, 3);
    else
        nchannels = 1;
    end
elseif ndims(images_raw) == 4
    nchannels = size(images_raw, 3);
else
    nchannels = 1;
end

channels = find(channel_mean_center);

h = waitbar(0, 'Applying image functions...');

for i=1:nimages
    waitbar(i/nimages, h);
    
    if dim == 3
        if nchannels > 1
            im_tmp = images(:,:,:,:,i);
            for i2=1:nstack
                for j=1:nchannels
                    im_tmp(:,:,j, i2) = apply_image_functions_oneimage(im_tmp(:,:,j, i2), channel_normalize(j), channel_blur(j), channel_weight(j));
                end
            end
            if ~isempty(channels)
                im_tmp = crop_zstack(im_tmp, channels, resize_image);
            end
            images(:,:,:,:,i) = im_tmp;
        else
            im_tmp = images(:,:,:,i);
            for i2=1:nstack
                % extract relevant channel; resize to number of pixels
                im_tmp(:,:, i2) = apply_image_functions_oneimage(im_tmp(:,:,i2), channel_normalize, channel_blur, channel_weight);
            end
            if channel_mean_center == 1
                im_tmp = crop_zstack(im_tmp, 0, resize_image);
            end
            images(:,:,:,i) = im_tmp;
        end
    else
        if nchannels > 1
            im_tmp = images(:,:,:,i);
            for j=1:nchannels
                im_tmp(:,:,j) = apply_image_functions_oneimage(im_tmp(:,:,j), channel_normalize(j), channel_blur(j), channel_weight(j));
            end
            if ~isempty(channels)
                im_tmp = crop_image(im_tmp, channels, resize_image);
            end
            images(:,:,:,i) = im_tmp;   
        else 
            im_tmp = images(:,:,i);
            im_tmp = apply_image_functions_oneimage(im_tmp, channel_normalize, channel_blur, channel_weight);
            if channel_mean_center == 1
                im_tmp = crop_image(im_tmp, 0, resize_image);
            end
            images(:,:,i) = im_tmp;
        end
    end
end

close(h);

function im2 = apply_image_functions_oneimage(im1, normalize_signal, signal_blur_radius, signal_scale)

if normalize_signal
    im1 = adapthisteq(im1);
end

npixels = size(im1, 1);
blur_size = round(signal_blur_radius*npixels);
if blur_size > 0
    filt = fspecial('disk', blur_size);
    im1 = imfilter(im1, filt, 'replicate');
end

im2 = immultiply(im1, signal_scale);
