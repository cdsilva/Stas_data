% read in dpERK pictures; align using vector dmaps

clear all
close all

res = '-r300';
fmt = '-djpeg';

cm_green = gray;
cm_green(:,1) = 0;
cm_green(:,3) = 0;

% % set parallel parameters
% run_in_parallel = true;
% nprocs = 2;
% 
% image_dir = '../membrane2/dpERK_staining';
% image_channel = 2;
% npixels = 100;
% max_deg = 15;
% 
% % load in membrane length data
% load('../membrane2/length_oct16.mat');
% mem_lengths = length(:,2);
% clear length;
% 
% %set image parameters
% %nimages = length(mem_lengths);
% nimages = 9;
% subplot_dim1 = floor(sqrt(nimages));
% subplot_dim2 = ceil(nimages / subplot_dim1);
% 
% [~, idx] = sort(mem_lengths, 'descend');
% mem_lengths = mem_lengths(idx(1:nimages));
% 
% % dimension of rotations
% rot_dim = 3;
% 
% % where to store fourier expansion coefficients
% ncoeff = sum(2*(0:max_deg)+1);
% F_l = zeros(nimages, ncoeff);
% 
% image_set = zeros(npixels, npixels, nimages);
% 
% % load in images
% figure;
% for i=1:nimages
%     im1 = imread(sprintf('%s/emb%02d.tif', image_dir, idx(i)));
%     %im1 = imread(sprintf('%s/emb%02d.tif', image_dir, i));
%     im1 = imresize(im1, [npixels npixels]);
%     im1 = im1(:,:,image_channel);
%     
%     subplot(subplot_dim1,subplot_dim2,i)
%     imshow(im1)
%     im1 = double(im1);
%     image_set(:,:,i) = im1;
%     
%     temp_F_l = zeros(1, ncoeff);
%     
%     for l=0:max_deg
%         start_idx = sum(2*(0:(l-1))+1)+1;
%         temp_F_l(start_idx:start_idx+2*l) = f_hat(l, im1);
%     end
%     
%     F_l(i,:) = temp_F_l;
% end
% 
% %% plot rotation of images
% [x, y, z] = sphere(2*npixels);
% theta = acos(z);
% phi = atan2(y, x);
% f = zeros(size(theta));
% 
% % test_ind = 1;
% % for l=0:max_deg
% %     start_idx = sum(2*(0:(l-1))+1)+1;
% %     sph_harm = Y_lm(l, reshape(theta, 1, []), reshape(phi, 1, []));
% %     f = f + reshape(sum(repmat(F_l(test_ind,start_idx:start_idx+2*l).', 1, size(sph_harm,2)).*conj(sph_harm), 1), size(f));
% % end
% % 
% % figure;
% % surf(x, y, z, real(f), 'linestyle','none')
% % colormap(gray)
% % 
% % 
% % [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, [0 pi/4 0]);
% % figure;
% % surf(x_rot, y_rot, z_rot, real(f), 'linestyle','none')
% % colormap(gray)
% % 
% % [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, [0 pi/2 0]);
% % figure;
% % surf(x_rot, y_rot, z_rot, real(f), 'linestyle','none')
% % colormap(gray)
% % 
% % [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, [0 3*pi/4 0]);
% % figure;
% % surf(x_rot, y_rot, z_rot, real(f), 'linestyle','none')
% % colormap(gray)
% % 
% % return
% 
% %% compute pairwise alignments
% 
% H_flat = zeros(nimages^2, rot_dim^2);
% 
% W_flat = zeros(nimages^2, 1);
% 
% opt_angles_flat = zeros(nimages^2, rot_dim);
% 
% if run_in_parallel
%     matlabpool('open',nprocs);
% end
% 
% siz = [nimages, nimages];
% 
% parfor ind=1:nimages^2
%     ind
%     
%     [i, j] = ind2sub(siz, ind);
%     
%     if i > j
%         tic
%         [F_l_aligned, A] = align_fhat(F_l(j,:), F_l(i,:), max_deg);
%         
%         opt_angles_flat(ind, :) = A;
%         
%         R = rot_matrix(A);
%         
%         H_flat(ind,:) = reshape(R, 1, []);
%         
%         W_flat(ind) = norm(F_l_aligned-F_l(i,:))^2;
%         toc
%     end
% end
% 
% 
% if run_in_parallel
%     matlabpool close;
% end
% 
% save('opt_alignments.mat','H_flat','W_flat','opt_angles_flat','F_l','mem_lengths','nimages','im_scale','max_deg');
% return;

load opt_alignments

%% reshape
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
[V, D] = eigs(H, 3);

% find optimal rotations
[u, s, v] = svd(V(1:3,1:3));
R0 = u*v';

image_set_aligned = zeros(size(image_set));

figure;
for i=1:nimages
    [u, s, v] = svd(V(3*i-2:3*i,1:3));
    R = u*v';
   
    im1 = image_set(:,:,i);

    [x2, y2, z2, f2] = project_full_image(im1);
    
%     if i == 1
%         A = get_angles([-1 0 0; 0 -1 0; 0 0 1]);
%     else
%         A = get_angles(R0*R');
%     end
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

figure;
plot(V(:,2),mem_lengths,'.')

[~, I] = sort(V(:,2));
figure;
for i=1:nimages
    
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set_aligned(:,:,I(i))))

end

%% vector DMAPS
eps = median(W(:));
W1 = exp(-W/eps);
W1 = diag(1./sum(W1)) * W1;
A = zeros(size(H));
for i=1:nimages
    for j=1:nimages
        A(3*i-2:3*i,3*j-2:3*j) = H(3*i-2:3*i,3*j-2:3*j) * W1(i,j);
    end
end

[V, D] = eigs(A, 5);

%% find optimal rotations

[u, s, v] = svd(V(1:3,1:3));
R0 = u*v';

image_set_aligned_vdm = zeros(size(image_set));

figure;
for i=1:nimages
    [u, s, v] = svd(V(3*i-2:3*i,1:3));
    R = u*v';
   
    im1 = image_set(:,:,i);

    [x2, y2, z2, f2] = project_full_image(im1);
    
%     if i == 1
%         A = get_angles([-1 0 0; 0 -1 0; 0 0 1]);
%     else
%         A = get_angles(R0*R');
%     end
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
n_embed = 4;
for embed_idx1=1:n_embed
    for embed_idx2 = 1:n_embed

        embed_coord = sum(reshape(V(:,embed_idx1), rot_dim, []) .* reshape(V(:,embed_idx2), rot_dim, []))';
        
        subplot(n_embed, n_embed, n_embed*(embed_idx1-1)+embed_idx2)
        plot(embed_coord,mem_lengths,'.')
        title(sprintf('i = %d, j = %d', embed_idx1, embed_idx2))
    end
end
    
    
plot(embed_coord,mem_lengths,'.')
%     end
% end

embed_idx1 = 1;
embed_idx2 = 1;

embed_coord = sum(reshape(V(:,embed_idx1), rot_dim, []) .* reshape(V(:,embed_idx2), rot_dim, []))';

figure;
plot(embed_coord,mem_lengths,'.')


[~, I] = sort(embed_coord);
figure;
for i=1:nimages
    
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(uint8(image_set_aligned(:,:,I(i))))

end


