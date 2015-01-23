function plot_images(images, image_dim, subplot_dim1, subplot_dim2)

nimages = size(images, ndims(images));
figure;
for i=1:nimages
    subplot(subplot_dim1, subplot_dim2, i)
    if image_dim == 3
        if ndims(images) == 4
            imshow(max_proj(images(:,:,:,i)))
        else
            imshow(max_proj(images(:,:,:,:,i)))
        end
    else
        if ndims(images) == 3
            imshow(images(:,:,i))
        else
            imshow(images(:,:,:,i))
        end
    end
end