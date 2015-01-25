function image2 = make_gray_nuclei(image)

image2 = image;

scale = 1.25;

image2(:,:,1) = immultiply(image2(:,:,1), scale);
for i=2:3
    image2(:,:,i) = imlincomb(1, image2(:,:,i), scale, image2(:,:,1));
end

image2 = circshift(image2, [0 0 -1]);

