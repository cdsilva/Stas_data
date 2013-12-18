function image_template = adjust_image_radially(image1)


image_template = image1;

[m, n] = size(image_template);

image_template = imadjust(image_template);

[X, Y] = meshgrid(1:m,1:n);
theta = atan2(Y-n/2, X-m/2);

delta = pi/20;
for i=-pi:delta:pi-delta
    data_IDX = find(theta <= i + delta & theta > i);
    image_template(data_IDX) = imadjust(image_template(data_IDX));
end