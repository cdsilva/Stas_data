
% set scattering transform parameters
filt_opt = struct();
filt_rot_opt = struct();
% oversampling is set to infty
% so that scattering of original and rotated
% image will be sampled at exactly the same points
scat_opt = struct();
scat_opt.oversampling = 10;

% compute scattering transform of first image
x = image_set_membrane(:,:,1);
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
sx_all = zeros(m, size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:m
    x = image_set_membrane(:,:,i);
    
    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

% dmaps

W = squareform(pdist(sx_all)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 20);

for i=3:20;
    figure;
    plot(V(:,2),V(:,i),'.')
end

for i=2:20
    corr(V(:,i),L(:,1))
end

idx = 8;

if corr(V(:,idx), L(:,1)) < 0
    V(:,idx) = -V(:,idx);
end

figure;
plot(L(:,1), V(:,idx),'.')
xlabel('membrane thickness')
ylabel(sprintf('$\\phi_%d$',idx),'interpreter','latex')
if print_figures
    print('DMAPS_membrane_scat_time_corr',fmt,res)
end

fprintf('Scattering transform (membrane) Spearman coeff: %2.4f \n', corr(L(:,1), V(:,idx), 'type','spearman'));


[~, I] = sort(V(:,idx));
if print_figures
    figure;
    for i=im_save_idx
        %subplot(subplot_dim1,subplot_dim2,i)

        imshow(uint8(image_set_membrane(:,:,I(i))), 'InitialMagnification', 'fit')
        %imshow(imadjust(uint8(image_set(:,:,I(i)))), 'InitialMagnification', 'fit')
        % make green colormap
        cm_green = gray;
        cm_green(:,1) = 0;
        cm_green(:,3) = 0;
        colormap(cm_green)
        axis off
        print(sprintf('membrane_scat_%d',i),fmt,res)
        clf
    end
end