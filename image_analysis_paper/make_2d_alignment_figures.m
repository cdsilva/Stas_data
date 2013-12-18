addpath('../membrane_pictures/synchron_so3')

load('opt_alignments.mat','H_flat','W_flat','nimages','rot_dim','image_set','subplot_dim1','subplot_dim2','siz','npixels')

%% reshape data
H = zeros(rot_dim*nimages);
W = zeros(nimages);

for i=1:nimages
    for j=1:i-1
        ind = sub2ind(siz, i, j);
        R = reshape(H_flat(ind,:), rot_dim, rot_dim);
        H((i-1)*rot_dim+1:i*rot_dim, (j-1)*rot_dim+1:j*rot_dim) = R;
        H((j-1)*rot_dim+1:j*rot_dim, (i-1)*rot_dim+1:i*rot_dim) = R';
        W(i,j) = W_flat(ind);
        W(j,i) = W_flat(ind);
    end
end
        
%% angular synchronization
R_opt = ang_synch(H, rot_dim);

R0 = R_opt(1:rot_dim, 1:rot_dim);

image_set_aligned = zeros(size(image_set));

figure;
for i=1:nimages
    R = R_opt(rot_dim*(i-1)+1:rot_dim*i, 1:rot_dim);
   
    im1 = image_set(:,:,i);

    [x2, y2, z2, f2] = project_full_image(im1);

    if i == 1
        f_rot = f2;
    else
        f_rot = rotate_fn2(f2, x2, y2, z2, R0*R');
    end
    
    im2 = extract_full_image(f_rot, [npixels npixels]);
    image_set_aligned(:,:,i) = im2;
    
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(im2))

end

W_new = zeros(nimages);

for i=1:nimages
    for j=1:i-1
        W_new(i,j) = sum(sum((image_set_aligned(:,:,i)-image_set_aligned(:,:,j)).^2));
        W_new(j,i) = W_new(i,j);
    end
end

eps = median(W_new(:));
[V, D] = dmaps(W_new, eps, 10);

if corr(V(:,2), L(:,1)) < 0
    V(:,2) = -V(:,2);
end

figure;
plot(L(:,1), V(:,2),'.')
xlabel('membrane thickness')
ylabel('\phi_2')
if print_figures
    print('angsynch_2d_time_corr',fmt, res)
end

[~, I] = sort(V(:,2));
if print_figures
    figure;
    for i=1:nimages  
        imshow(uint8(image_set_aligned(:,:,I(i))), 'InitialMagnification', 'fit')
        % make green colormap
        cm_green = gray;
        cm_green(:,1) = 0;
        cm_green(:,3) = 0;
        colormap(cm_green)
        axis off
        print(sprintf('dpERK_angsynch_%d',i),fmt,res)
        clf
    end
end

%% vector DMAPS
eps = median(W(:));
neigs = 6;
[R_opt, embed_coord, embed_idx, D] = vdm(H, W, eps, neigs);

R0 = R_opt(1:rot_dim, :);

image_set_aligned_vdm = zeros(size(image_set));

figure;
for i=1:nimages
    R = R_opt(rot_dim*(i-1)+1:i*rot_dim, :);
   
    im1 = image_set(:,:,i);

    [x2, y2, z2, f2] = project_full_image(im1);
    
    if i == 1
        f_rot = f2;
    else
        f_rot = rotate_fn2(f2, x2, y2, z2, R0*R');
    end
    
    im2 = extract_full_image(f_rot, [npixels npixels]);
    image_set_aligned_vdm(:,:,i) = im2;
    
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(im2))

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
    
coord_idx = 1;
if corr(embed_coord(:,coord_idx), L(:,1)) < 0
    embed_coord(:,coord_idx) = -embed_coord(:,coord_idx);
end

figure;
plot(L(:,1),embed_coord(:,coord_idx),'.')
xlabel('membrane thickness')
ylabel(sprintf('\\langle \\phi_%d, \\phi_%d \\rangle', embed_idx(1, coord_idx), embed_idx(2, coord_idx)))
if print_figures
    print('vdm_2d_time_corr',fmt, res)
end

[~, I] = sort(embed_coord(:,coord_idx));
if print_figures
    figure;
    for i=1:nimages  
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

rmpath('../membrane_pictures/synchron_so3')