function im_new = shift_image(im1, alpha, beta, gamma, angle_proj)

im_new = im1;

im_new = rotate_image(im_new, alpha);
im_new = translate_image(im_new, beta, gamma, angle_proj);