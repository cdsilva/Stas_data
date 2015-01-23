function images2 = order_all_images(images, embed_coord)

[~, I] = sort(embed_coord);

if ndims(images) == 3
    images2 = images(:,:,I);
elseif ndims(images) == 4
    images2 = images(:,:,:,I);
elseif ndims(images) == 5
    images2 = images(:,:,:,:,I);
end