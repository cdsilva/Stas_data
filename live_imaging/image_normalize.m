function image2 = image_normalize(image, rot_angle)

npixels = 100;
crop_threshold = 1e3;

image_resize = @(image) imresize(image, [npixels npixels]);
image_equalize = @(image, rot_angle) crop_image(adapthisteq(imrotate(image, rot_angle, 'crop'), 'cliplimit',0.005), crop_threshold, 0);
image_symmetrize = @(image) imlincomb(0.5, image, 0.5, fliplr(image));

image2 = adapthisteq(image_symmetrize(image_equalize(image_resize(image), rot_angle)));