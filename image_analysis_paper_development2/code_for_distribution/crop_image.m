function image2 = crop_image(image, channel, resize_image)
% image2 = crop_image(image, channel)
% crop the black border from an image
% image is the image to be cropped
% channel is the channel to use for the border detection for an RGB image
% (if the image is grayscale, the value of channel is ignored)
% image2 is the cropped image, with constant 10% border added around the edges

% threshold used to detect edges
edge_tol = 0;
% type of edge filer
detection_type = 'canny';

npixels = size(image, 1);

% compute edges in image
if ndims(image) == 3
    im_tmp = image(:,:,channel);
    if length(channel) > 1
        im_tmp = uint8(mean(double(im_tmp), 3));
    end
else
    im_tmp = image;
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
    if ndims(image) == 3
        image2 = padarray(image, [nbuffer nbuffer 0]);
        image2 = imresize(image2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),:), [npixels npixels]);
    else
        image2 = padarray(image, [nbuffer nbuffer]);
        image2 = imresize(image2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer)), [npixels npixels]);
    end
else
    RI = imref2d([npixels npixels],[1 npixels],[1 npixels]);
    tform = affine2d([1 0 0; 0 1 0; (1+npixels)/2-(idx1+idx2)/2 (1+npixels)/2-(idx3+idx4)/2 1]);
    image2 = imwarp(image, RI, tform, 'outputview', RI);
end
