function save_images(images, dim, image_dir, image_name, image_ext, stack_name)

nimages = size(images, ndims(images));

if ~exist(image_dir, 'dir')
    mkdir(image_dir);
end

h = waitbar(0, 'Reading images...');
if dim == 2
    if ndims(images) == 3
        for i=1:nimages
            waitbar(i/nimages, h);
            filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
            write_single_image(images(:,:,i), filename)
        end
    else
        for i=1:nimages
            waitbar(i/nimages, h);
            filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
            write_single_image(images(:,:,:,i), filename)
        end
    end
    
elseif dim == 3
    nstack = size(images, ndims(images)-1);
    if ndims(images) == 4
        for i=1:nimages
            waitbar(i/nimages, h);
            if ~exist(sprintf('%s/%s%02d', image_dir, image_name, i), 'dir')
                mkdir(sprintf('%s/%s%02d', image_dir, image_name, i));
            end
            for j=1:nstack
                filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
                write_single_image(images(:,:,j,i), filename);
            end
        end
    else
        for i=1:nimages
            waitbar(i/nimages, h);
            if ~exist(sprintf('%s/%s%02d', image_dir, image_name, i), 'dir')
                mkdir(sprintf('%s/%s%02d', image_dir, image_name, i));
            end
            for j=1:nstack
                filename = sprintf('%s/%s%02d/%s%02d.%s', image_dir, image_name, i, stack_name, j, image_ext);
                write_single_image(images(:,:,:,j,i), filename);
            end
        end
    end
end
close(h);

function write_single_image(im1, filename)

if exist(filename, 'file')
    choice = questdlg(sprintf('%s already exists. Would you like to overwrite the image?', filename), ...
        'File already exists', ...
        'Yes',...
        'No',...
        'No');
    % Handle response
    switch choice
        case 'Yes'
            imwrite(im1,filename);
        case 'No'
            return
    end
else
    imwrite(im1,filename);    
end