addpath('../membrane_pictures/synchron_so3_nosphharm')

%% parameters

% amont of border or buffer to add around the images
buffer_size = 20;

% maximum shift in pixels
shift_max = 20;

%size of "portion" of sphere on which to project
angle_proj = pi/4;

% total number of pixels
npixels2 = npixels+2*buffer_size;

% dimension of rotations
dim = 3;

% store images in 3d array
image_set_buffered = zeros(npixels2, npixels2, m);

%% load in images

for i=1:m
    %store image
    im1 = image_set(:,:,i);
    image_set_buffered(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels,i) = im1;
end

%% compute pairwise alignments
tic
xmax = 0.1;
ymax = 0.1;
[R, W, angles] = align_data_nosph(image_set_buffered, xmax, ymax, shift_max);

        
%% angular synchronization
% dim = 3;
% 
% R_opt = ang_synch(R, dim);
% 
% image_set_aligned = zeros(size(image_set));
% 
% figure;
% for i=1:m
%     
%     image_tmp = image_set_buffered(:,:,i);
%     
%     R_all = R_opt(1:3, 1:3)'*R_opt(3*i-2:3*i,:);
%     
%     [alpha, beta, gamma] = calc_angles(R_all);    
%     image_tmp = shift_image(image_tmp, alpha, beta, gamma, angle_proj);
%     image_set_aligned(:,:,i) = image_tmp(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels);
%     
%     subplot(subplot_dim1, subplot_dim2, i);
%     imshow(uint8(image_set_aligned(:,:,i)),'InitialMagnification', 'fit')
% 
% end
% 
% W_new = zeros(m);
% 
% for i=1:m
%     for j=1:i-1
%         W_new(i,j) = sum(sum((image_set_aligned(:,:,i)-image_set_aligned(:,:,j)).^2));
%         W_new(j,i) = W_new(i,j);
%     end
% end
% 
% eps = median(W_new(:));
% [V, D] = dmaps(W_new, eps, 10);
% 
% if corr(V(:,2), L(:,1)) < 0
%     V(:,2) = -V(:,2);
% end
% 
% figure;
% plot(L(:,1), V(:,2),'.')
% xlabel('membrane thickness')
% ylabel('$\phi_2$','interpreter','latex')
% if print_figures
%     print('angsynch_2d_time_corr',fmt, res)
% end
% 
% [~, I] = sort(V(:,2));
% if print_figures
%     figure;
%     for i=im_save_idx
%         imshow(uint8(image_set_aligned(:,:,I(i))), 'InitialMagnification', 'fit')
%         % make green colormap
%         cm_green = gray;
%         cm_green(:,1) = 0;
%         cm_green(:,3) = 0;
%         colormap(cm_green)
%         axis off
%         print(sprintf('dpERK_angsynch_%d',i),fmt,res)
%         clf
%     end
% end

%% vector DMAPS

eps = median(W(:));
neigs = 5;
[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
runtime = toc

image_set_aligned_vdm = zeros(size(image_set));

figure;
for i=1:m
    image_tmp = image_set_buffered(:,:,i);
    
    R_all = (R_opt(1:3, 1:3)'*R_opt(3*i-2:3*i,:));
    
    [theta, xdisp, ydisp] = calc_angles(R_all, xmax, ymax, npixels2);
    theta
    xdisp
    ydisp
    
    image_tmp = shift_image2(image_tmp, theta, xdisp, ydisp);
    image_set_aligned_vdm(:,:,i) = image_tmp(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels);
    
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_set_aligned_vdm(:,:,i)),'InitialMagnification', 'fit')

end

return

figure;
for i=1:m
    image_tmp = image_set_buffered(:,:,i);
    
    R_all = R(3*i-2:3*i,1:3);
    
    [theta, xdisp, ydisp] = calc_angles(R_all, xmax, ymax, npixels2);
    theta
    
    image_tmp = shift_image2(image_tmp, theta, xdisp, ydisp);
    image_set_aligned_vdm(:,:,i) = image_tmp(buffer_size+1:buffer_size+npixels,buffer_size+1:buffer_size+npixels);
    
    subplot(subplot_dim1, subplot_dim2, i);
    imshow(uint8(image_set_aligned_vdm(:,:,i)),'InitialMagnification', 'fit')

end


figure;
n_embed = size(embed_coord, 2);
splot_1 = ceil(sqrt(n_embed));
splot_2 = ceil(n_embed / splot_1);
for i=1:n_embed   
    subplot(splot_1, splot_2, i)
    plot(embed_coord(:,i),L(:,1),'.')
    title(sprintf('i = %d, j = %d', embed_idx(1,i), embed_idx(2,i)))
end
    
[~, coord_idx] = max(abs(corr(embed_coord, L(:,1))));
if corr(embed_coord(:,coord_idx), L(:,1)) < 0
    embed_coord(:,coord_idx) = -embed_coord(:,coord_idx);
end

[~, I] = sort(embed_coord(:,coord_idx));

% figure;
% plot(L(:,1),embed_coord(:,coord_idx),'.')
% hold on
% plot(L(I(im_save_idx),1), embed_coord(I(im_save_idx),coord_idx), 'o')
% xlabel('membrane thickness')
% ylabel(sprintf('$\\langle \\phi_%d, \\phi_%d \\rangle$', embed_idx(1, coord_idx), embed_idx(2, coord_idx)),'interpreter','latex')
% if print_figures
%     print('vdm_2d_time_corr',fmt, res)
% end

fprintf('VDM (2-D) Spearman coeff: %2.4f \n', corr(L(:,1),embed_coord(:,coord_idx), 'type','spearman'));



if print_figures
    figure;
    for i=im_save_idx 
        set(gcf, 'paperposition',[0 0 8 8])
        imshow(uint8(image_set_aligned_vdm(:,:,I(i))), 'InitialMagnification', 'fit')
        % make green colormap
        cm_green = gray;
        cm_green(:,1) = 0;
        cm_green(:,3) = 0;
        colormap(cm_green)
        axis off
        print(sprintf('dpERK_vdm_%d',i),fmt,res)
        clf
    end
end

figure;
for i=1:m
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(uint8(image_set_aligned_vdm(:,:,I(i))), 'InitialMagnification', 'fit')
    %set(gca,'position',[0 0 1 1],'units','normalized')
    % make green colormap
    cm_green = gray;
    cm_green(:,1) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
end

%save('2d_alignment_figures.mat');

%%
rmpath('../membrane_pictures/synchron_so3_nosphharm')