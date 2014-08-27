function image2 = image_fn_shuffledmovie(image, npixels)

filt = fspecial('disk');

image2 = imresize(image, [npixels npixels]);

channel = 0;
image2 = crop_image(image2, channel);

image2 = padarray(image2, [10 10]);
image2 = imresize(image2, [npixels npixels]);

image2 = adapthisteq(image2);
% tmp_vec = round(linspace(1, npixels, 3));
% for i=1:length(tmp_vec)-1
%     for j=1:length(tmp_vec)-1
%         image2(tmp_vec(i):tmp_vec(i+1), tmp_vec(j):tmp_vec(j+1)) = imadjust(image2(tmp_vec(i):tmp_vec(i+1), tmp_vec(j):tmp_vec(j+1)));
%     end
% end

image2 = imfilter(image2, filt, 'replicate');