function image2 = image_fn(image, npixels)

filt = fspecial('disk');

image2 = imresize(image, [npixels npixels]);

channel = 0;
image2 = crop_image(image2, channel);

image2 = padarray(image2, [10 10]);
image2 = imresize(image2, [npixels npixels]);

image2 = adapthisteq(image2);

image2 = imfilter(image2, filt, 'replicate');