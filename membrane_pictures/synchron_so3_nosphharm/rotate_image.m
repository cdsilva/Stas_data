function im_new = rotate_image(im1, alpha)

im_new = im1;
im_new = imrotate(im_new, alpha*(180/pi), 'crop');
