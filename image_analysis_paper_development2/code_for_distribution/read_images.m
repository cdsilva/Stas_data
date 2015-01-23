function [images, nchannels] = read_images(image_dir, image_name, image_ext, stack_name, nimages, nstack, npixels, dim)

try
    h = waitbar(0, 'Reading images...');
    if dim == 2
        i = 1;
        filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
        im_tmp = imread(filename);
        if ndims(im_tmp) == 2
            nchannels = 1;
            
            images = zeros(npixels, npixels, nimages, 'uint8');
            for i=1:nimages
                waitbar(i/nimages, h);
                filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
                im_tmp = imread(filename);
                im_tmp = imresize(im_tmp, [npixels npixels]);
                
                images(:, :, i) = im_tmp;
            end
        else
            nchannels = size(im_tmp, 3);
            
            images = zeros(npixels, npixels, nchannels, nimages, 'uint8');
            for i=1:nimages
                waitbar(i/nimages, h);
                filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
                im_tmp = imread(filename);
                im_tmp = imresize(im_tmp, [npixels npixels]);
                
                images(:, :, :, i) = im_tmp;
            end
        end
    elseif dim == 3
        
        i = 1;
        j = 1;
        filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
        im_tmp = imread(filename);
        if ndims(im_tmp) == 2
            nchannels = 1;
            
            images = zeros(npixels, npixels, nstack, nimages, 'uint8');
            for i=1:nimages
                waitbar(i/nimages, h);
                for j=1:nstack
                    filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
                    im_tmp = imread(filename);
                    im_tmp = imresize(im_tmp, [npixels npixels]);
                    
                    images(:, :, j, i) = im_tmp;
                end
            end
        else
            nchannels = size(im_tmp, 3);
            
            images = zeros(npixels, npixels, nchannels, nstack, nimages, 'uint8');
            for i=1:nimages
                waitbar(i/nimages, h);
                for j=1:nstack
                    filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
                    
                    im_tmp = imread(filename);
                    im_tmp = imresize(im_tmp, [npixels npixels]);
                    
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