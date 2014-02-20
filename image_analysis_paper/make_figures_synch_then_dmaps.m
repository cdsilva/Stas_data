clear all
close all

load dpERK_aligned_small_all

R_opt = ang_synch(R, dim);

image_set_aligned = zeros(size(image_set));
for i=1:m
    im1 = rotate_image(R_opt(dim*(i-1)+1:dim*i,:), uint8(image_set(:,:,i)), angle_proj);
    image_set_aligned(:,:,i) = double(im1);
end

W = zeros(m);
for i=1:m
    for j=1:i-1
        W(i,j) = sum(sum((image_set_aligned(:,:,i)-image_set_aligned(:,:,j)).^2));
    end
end

W = W+W';

[V, D] = dmaps(W, eps, 10);

if corr(V(:,2),mem_lengths) < 0
    V(:,2) = - V(:,2);
end

[~, I] = sort(V(:,2));
figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_set_aligned(:,:,I(i))));
end
print('synch_then_dmaps_images','-r300','-djpeg');


ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_vdm = compute_ranks(V(:,2));


figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('rank from membrane lengths')
ylabel('rank from dmaps')
print('synch_then_dmaps_corr','-r300','-djpeg');

