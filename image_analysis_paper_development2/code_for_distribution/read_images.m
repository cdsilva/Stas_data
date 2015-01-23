function [images, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim)
if dim == 2
    
    i = 1;
    im_tmp = imread(sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext));
    if ndims(im_tmp) == 2
        nchannels = 1;
        
        images = zeros(npixels, npixels, nimages, 'uint8');
        for i=1:nimages
            
            im_tmp = imread(sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext));
            im_tmp = imresize(im_tmp, [npixels npixels]);
            
            images(:, :, i) = im_tmp;
        end
    else
        nchannels = size(im_tmp, 3);
        
        images = zeros(npixels, npixels, nchannels, nimages, 'uint8');
        for i=1:nimages
            
            im_tmp = imread(sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext));
            im_tmp = imresize(im_tmp, [npixels npixels]);
            
            images(:, :, :, i) = im_tmp;
        end
    end
elseif dim == 3
    
    i = 1;
    j = 1;
    im_tmp = imread(sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext));
    if ndims(im_tmp) == 2
        nchannels = 1;
        
        images = zeros(npixels, npixels, nstack, nimages, 'uint8');
        for i=1:nimages
            for j=1:nstack
                
                im_tmp = imread(sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext));
                im_tmp = imresize(im_tmp, [npixels npixels]);
                
                images(:, :, j, i) = im_tmp;
            end
        end
    else
        nchannels = size(im_tmp, 3);
        
        images = zeros(npixels, npixels, nchannels, nstack, nimages, 'uint8');
        for i=1:nimages
            for j=1:nstack
                
                im_tmp = imread(sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext));
                im_tmp = imresize(im_tmp, [npixels npixels]);
                
                images(:, :, :, j, i) = im_tmp;
            end
        end
    end
else
    h = msgbox('ERROR: invalid image dimesion');
    return
end