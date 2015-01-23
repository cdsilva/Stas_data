function image_stack2 = crop_zstack(image_stack, channel)

im_tmp = uint8(mean(double(image_stack), ndims(image_stack)));
if ndims(im_tmp) == 3
        im_tmp = im_tmp(:,:,channel);
%     im_tmp = image;
%     for i=1:size(image, 3)
%         im_tmp(:,:,i) = imadjust(im_tmp(:,:,i));
%     end
    if length(channel) > 1
        im_tmp = uint8(mean(double(im_tmp), 3));
    end
end
im_tmp = imadjust(im_tmp);

npixels = size(im_tmp, 1);

bw = activecontour(im_tmp, im2bw(im_tmp, 0.1));
bw = imfill(bw, 'holes');
% STATS = regionprops(bw, 'Centroid');

[y x] = find( bw );

RI = imref2d([npixels npixels],[1 npixels],[1 npixels]);
% tform = affine2d([1 0 0; 0 1 0; (1+npixels)/2-STATS.Centroid(1) (1+npixels)/2-STATS.Centroid(2) 1]);
tform = affine2d([1 0 0; 0 1 0; (1+npixels)/2-mean(x) (1+npixels)/2-mean(y) 1]);

image_stack2 = imwarp(image_stack, RI, tform, 'outputview', RI);
