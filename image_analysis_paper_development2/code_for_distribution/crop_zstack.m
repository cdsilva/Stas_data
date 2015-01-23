function image_stack2 = crop_zstack(image_stack, channel, resize_image)

% threshold used to detect edges
edge_tol = 0;
% type of edge filer
detection_type = 'canny';

npixels = size(image_stack, 1);

im_tmp = uint8(mean(double(image_stack), ndims(image_stack)));
if ndims(im_tmp) == 3
    im_tmp = im_tmp(:,:,channel);
    if length(channel) > 1
        im_tmp = uint8(mean(double(im_tmp), 3));
    end
end
im_tmp = imadjust(im_tmp);

im_tmp = edge(im_tmp, detection_type);

% find indices where edges start and end
row_sum = sum(im_tmp);
idx1 = find(row_sum > edge_tol, 1, 'first');
idx2 = find(row_sum > edge_tol, 1, 'last');

col_sum = sum(im_tmp, 2);
idx3 = find(col_sum > edge_tol, 1, 'first');
idx4 = find(col_sum > edge_tol, 1, 'last');

% crop image; add border of constant width
if resize_image
    nbuffer = round(0.1*npixels);
    if ndims(image_stack) == 4
        image_stack2 = padarray(image_stack, [nbuffer nbuffer 0 0]);
        image_stack2 = imresize(image_stack2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),:,:), [npixels npixels]);
    else
        image_stack2 = padarray(image_stack, [nbuffer nbuffer 0]);
        image_stack2 = imresize(image_stack2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),:), [npixels npixels]);
    end
else
    RI = imref2d([npixels npixels],[1 npixels],[1 npixels]);
    tform = affine2d([1 0 0; 0 1 0; (1+npixels)/2-(idx1+idx2)/2 (1+npixels)/2-(idx3+idx4)/2 1]);
    image_stack2 = imwarp(image_stack, RI, tform, 'outputview', RI);
end