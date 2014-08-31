function image2 = crop_image(image, channel)
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
    im_tmp = edge(image(:, :, channel), detection_type);
else
    im_tmp = edge(image, detection_type);
end

% find indices where edges start and end
row_sum = sum(im_tmp);
idx1 = find(row_sum > edge_tol, 1, 'first');
idx2 = find(row_sum > edge_tol, 1, 'last');

col_sum = sum(im_tmp, 2);
idx3 = find(col_sum > edge_tol, 1, 'first');
idx4 = find(col_sum > edge_tol, 1, 'last');

% crop image; add border of constant width
if ndims(image) == 3
    image2 = imresize(image(idx3:idx4, idx1:idx2,:), [npixels npixels]);
    image2 = padarray(image2, [round(0.1*npixels) round(0.1*npixels) 0]);
    image2 = imresize(image2, [npixels npixels]);
elseif ndims(image) == 2
    image2 = imresize(image(idx3:idx4, idx1:idx2), [npixels npixels]);    
    image2 = padarray(image2, [round(0.1*npixels) round(0.1*npixels)]);
    image2 = imresize(image2, [npixels npixels]);
end