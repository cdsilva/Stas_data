function [images, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim)
%READ_IMAGES reads in images stored in specified directory and
%subdirectories
% [images, nchannels] = READ_IMAGES(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim)
% reads in the images stored in image_dir
% image_name gives the image prefix (for 2D images) or folder prefix (for
% 3D z-stacks) for the images
% image_ext is the extension of the images (tif, jpg, etc.)
% stack_name is the image prefix for each of the images in a z-stack; if
% there are only 2D images, then stack_name is ignored
% nimages is the number of 2D images, or the number of 3D z-stacks
% nstack is the number of images in a single z-stack; if
% there are only 2D images, then stack_name is ignored
% npixels is the number of pixels for reading in the images (this does not
% have to be the same as the actual resolutoin of the images)
% dim is the image dimension: dim=2 for standard 2D images, dim=3 for
% z-stacks
%
% images is the array of images
% nchannels is the number of channels in the images

try
    h = waitbar(0, 'Reading images...');
    
    if dim == 2
        i = 1;
        filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
    else
        i = 1;
        j = 1;
        filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
    end
    im_tmp = imread(filename);
    if ndims(im_tmp) == 2
        nchannels = 1;
    else
        nchannels = size(im_tmp, 3);
    end
    
    if dim == 2
        if nchannels == 1
            images = zeros(npixels, npixels, nimages, 'uint8');
        else
            images = zeros(npixels, npixels, nchannels, nimages, 'uint8');
        end
    else
        if nchannels == 1
            images = zeros(npixels, npixels, nstack, nimages, 'uint8');
        else
            images = zeros(npixels, npixels, nchannels, nstack, nimages, 'uint8');
        end
    end
    
    for i=1:nimages
        waitbar(i/nimages, h);
        if dim == 2
            filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
            im_tmp = imread(filename);
            im_tmp = imresize(im_tmp, [npixels npixels]);
            
            if nchannels == 1
                images(:, :, i) = im_tmp;
            else
                images(:, :, :, i) = im_tmp;
            end
        else
            for j=1:nstack
                filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
                im_tmp = imread(filename);
                im_tmp = imresize(im_tmp, [npixels npixels]);
                
                if nchannels == 1
                    images(:, :, j, i) = im_tmp;
                else
                    images(:, :, :, j, i) = im_tmp;
                end
            end
        end
    end
    close(h);
    
catch ME
    if strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
        close(h);
        msgbox(sprintf('%s does not exist', filename));
        rethrow(ME);
    end
end