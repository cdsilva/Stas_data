function im_new = translate_image(im1, beta, gamma, angle_proj)

[n, n1] = size(im1);
if n ~= n1
    disp('Image is not square');
    return;
end

im_new = im1;

im_new = circshift(im_new, round([beta gamma] * n/angle_proj));