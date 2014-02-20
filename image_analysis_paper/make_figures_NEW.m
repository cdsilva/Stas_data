clear all
close all

res = '-r300';
fmt = '-djpeg';

im_save_dir = 'paper_figures';

%% dpERK
load dpERK_aligned_all
file_init = 'dpERK';

m = size(image_set_aligned2, 3);
dim1 = ceil(sqrt(m));
dim2 = ceil(m/dim1);

figure;
set(gcf, 'paperposition',[0 0 8 8])
cm_red = gray;
cm_red(:,2) = 0;
cm_red(:,3) = 0;
for i=1:m
    imshow(image_set_aligned2(:,:,i))
    set(gca,'position',[0 0 1 1],'units','normalized')
    colormap(cm_red);
    
    print(sprintf('%s/%s_%d', im_save_dir,file_init, i), fmt, res);
    clf
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    subplot(dim1, dim2, i)
    imshow(image_set_aligned2(:,:,i))
    colormap(cm_red);
end
print(sprintf('%s/%s_array', im_save_dir,file_init), fmt, res);

figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('ranking from ordering by membrane length', 'fontsize', 20)
ylabel('ranking from vdm using dpERK', 'fontsize', 20)
print(sprintf('%s/%s_rank_corr', im_save_dir, file_init), fmt, res);
disp(corr(ranks_from_membranes, ranks_from_vdm))

load('../membrane_pictures/large_dataset/dpERK04_big.mat');
dpERK = dpERK(ind, :);
figure;
imagesc(dpERK(I,:))
print(sprintf('%s/%s_onedimensional', im_save_dir, file_init), fmt, res);


%% Dl
load Dl_aligned_all
file_init = 'Dl';

m = size(image_set_aligned2, 3);
dim1 = ceil(sqrt(m));
dim2 = ceil(m/dim1);

figure;
set(gcf, 'paperposition',[0 0 8 8])
cm_green = gray;
cm_green(:,1) = 0;
cm_green(:,3) = 0;
for i=1:m
    imshow(image_set_aligned2(:,:,i))
    set(gca,'position',[0 0 1 1],'units','normalized')
    colormap(cm_green);
    
    print(sprintf('%s/%s_%d', im_save_dir,file_init, i), fmt, res);
    clf
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
for i=1:m
    subplot(dim1, dim2, i)
    imshow(image_set_aligned2(:,:,i))
    colormap(cm_green);
    
end
print(sprintf('%s/%s_array', im_save_dir,file_init), fmt, res);

figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('ranking from ordering by membrane length', 'fontsize', 20)
ylabel('ranking from vdm using Dl', 'fontsize', 20)
print(sprintf('%s/%s_rank_corr', im_save_dir, file_init), fmt, res);
disp(corr(ranks_from_membranes, ranks_from_vdm))

%% dpERK
load dpERK_aligned_small_all
file_init = 'dpERK_small';

m = size(image_set_aligned2, 3);
dim1 = ceil(sqrt(m));
dim2 = ceil(m/dim1);

figure;
set(gcf, 'paperposition',[0 0 8 8])
cm_red = gray;
cm_red(:,2) = 0;
cm_red(:,3) = 0;
for i=1:m
    subplot(dim1, dim2, i)
    imshow(image_set_aligned2(:,:,i))
    colormap(cm_red);
end
print(sprintf('%s/%s_array', im_save_dir,file_init), fmt, res);

figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('ranking from ordering by membrane length', 'fontsize', 20)
ylabel('ranking from vdm using dpERK', 'fontsize', 20)
print(sprintf('%s/%s_rank_corr', im_save_dir, file_init), fmt, res);
disp(corr(ranks_from_membranes, ranks_from_vdm))

load('../membrane_lengths/oct16.mat');
figure;
imagesc(dpERK_raw(I,:))
print(sprintf('%s/%s_onedimensional', im_save_dir, file_init), fmt, res);

%% Dl
load Dl_aligned_small_all
file_init = 'Dl_small';

m = size(image_set_aligned2, 3);
dim1 = ceil(sqrt(m));
dim2 = ceil(m/dim1);

figure;
set(gcf, 'paperposition',[0 0 8 8])
cm_green = gray;
cm_green(:,1) = 0;
cm_green(:,3) = 0;
for i=1:m
    subplot(dim1, dim2, i)
    imshow(image_set_aligned2(:,:,i))
    colormap(cm_green);
    
end
print(sprintf('%s/%s_array', im_save_dir,file_init), fmt, res);

figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('ranking from ordering by membrane length', 'fontsize', 20)
ylabel('ranking from vdm using Dl', 'fontsize', 20)
print(sprintf('%s/%s_rank_corr', im_save_dir, file_init), fmt, res);

disp(corr(ranks_from_membranes, ranks_from_vdm))

