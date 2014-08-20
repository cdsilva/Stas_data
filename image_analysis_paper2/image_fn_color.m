function image2 = image_fn_color(image, channel, npixels)

filt = fspecial('disk', 5);

image2 = imresize(image, [npixels npixels]);

image2 = crop_image(image2, channel);
image2 = padarray(image2, [10 10 0]);
image2 = imresize(image2, [npixels npixels]);

image2(:, :, channel) = adapthisteq(image2(:, :, channel));
image2(:, :, channel) = imfilter(image2(:, :, channel), filt, 'replicate');
% image2(:, :, channel) = imadjust(image2(:, :, channel));
image2(:, :, channel) = immultiply(image2(:, :, channel), 0.5);
