function image2 = image_normalize(image, rot_angle)


npixels = 100;

image2 = imresize(image, [npixels npixels]);
image2 = imrotate(image2, rot_angle, 'crop');
image2 = adapthisteq(image2, 'cliplimit',0.005);
image2 = crop_image(image2, 0);
image2 = imlincomb(0.5, image2, 0.5, fliplr(image2));
image2 = adapthisteq(image2);