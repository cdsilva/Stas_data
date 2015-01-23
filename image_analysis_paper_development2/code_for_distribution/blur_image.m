function image2 = blur_image(image, filter_width)

filter_size = round(filter_width*size(image, 1));
if filter_size > 0
    filt = fspecial('disk', filter_size);
    image2 = imfilter(image, filt, 'replicate');
end