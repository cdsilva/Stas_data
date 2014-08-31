function image2 = normalize_nuclear_signal(image)
% normalize_nuclear_signal(image)
% normalize the image which contains the nuclear signal, using adaptive
% histogram equalization + disk filter

npixels = size(image, 1);

filt = fspecial('disk', round(0.05*npixels));

image2 = adapthisteq(image);

image2 = imfilter(image2, filt, 'replicate');