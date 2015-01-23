function plot_images(images, image_dim)

nimages = size(images, ndims(images));

subplot_dim1 = round(sqrt(nimages));
subplot_dim2 = ceil(nimages/subplot_dim1);

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

function image1 = max_proj(zstack)

image1 = max(zstack, [], ndims(zstack));