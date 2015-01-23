function save_images(images, dim, image_dir, image_name, image_ext, stack_name)

nimages = size(images, ndims(images));

if ~exist(image_dir, 'dir')
    mkdir(image_dir);
end

h = waitbar(0, 'Saving images...');
if dim == 2
    if ndims(images) == 3
        for i=1:nimages
            waitbar(i/nimages, h);
            filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
            check_val = write_single_image(images(:,:,i), filename);
            if check_val == -1
                close(h);
                return
            end
        end
    else
        for i=1:nimages
            waitbar(i/nimages, h);
            filename = sprintf('%s/%s%02d.%s', image_dir, image_name, i, image_ext);
            check_val = write_single_image(images(:,:,:,i), filename);
            if check_val == -1
                close(h);
                return
            end
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
                check_val = write_single_image(images(:,:,j,i), filename);
                if check_val == -1
                    close(h);
                    return
                end
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
                check_val = write_single_image(images(:,:,:,j,i), filename);
                if check_val == -1
                    close(h);
                    return
                end
            end
        end
    end
end
close(h);

function check_val = write_single_image(im1, filename)

if exist(filename, 'file')
    choice = questdlg(sprintf('%s already exists. Would you like to overwrite the image?', filename), ...
        'File already exists', ...
        'Yes',...
        'No',...
        'Cancel writing files',...
        'No');
    % Handle response
    switch choice
        case 'Yes'
            check_val = 0;
            imwrite(im1,filename);
        case 'No'
            check_val = 0;
            return
        case 'Cancel writing files'
            check_val = -1;
            return
    end
else
    check_val = 0;
    imwrite(im1,filename);    
end