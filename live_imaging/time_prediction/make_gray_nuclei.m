function image2 = make_gray_nuclei(image)

scale = 0.5;

image2 = image;
image2(:,:,1) = 0;
image2 = imlincomb(scale, repmat(image(:,:,1), [1 1 3]), 1, image2);

image2 = circshift(image2, [0 0 -1]);

