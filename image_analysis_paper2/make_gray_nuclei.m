function image2 = make_gray_nuclei(image)

image2 = image;
for i=2:3
    image2(:,:,i) = imlincomb(1, image2(:,:,i), 1, image2(:,:,1));
end

image2 = circshift(image2, [0 0 -1]);

