function image2 = crop_image2(image)

edge_tol = 2;
detection_type = 'canny'; 

npixels = size(image, 1);

if ndims(image) == 3
    im_tmp = edge(imlincomb(1, image(:, :, 1), 1, image(:, :, 2), 1, image(:, :, 3)), detection_type);
else
    im_tmp = edge(image, detection_type);
end

row_sum = sum(im_tmp);
idx1 = find(row_sum > edge_tol, 1, 'first');
idx2 = find(row_sum > edge_tol, 1, 'last');

col_sum = sum(im_tmp, 2);
idx3 = find(col_sum > edge_tol, 1, 'first');
idx4 = find(col_sum > edge_tol, 1, 'last');

if ndims(image) == 3
    image2 = imresize(image(idx3:idx4, idx1:idx2,:), [npixels npixels]);
elseif ndims(image) == 2
    image2 = imresize(image(idx3:idx4, idx1:idx2), [npixels npixels]);    
end