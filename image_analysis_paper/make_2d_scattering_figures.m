%addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
addpath '../../../MATLAB/scatnet-0.2';
addpath_scatnet

%% compute scattering coefficients
% set scattering transform parameters
filt_opt = struct();
filt_rot_opt = struct();
% oversampling is set to infty
% so that scattering of original and rotated
% image will be sampled at exactly the same points
scat_opt = struct();
scat_opt.oversampling = 10;

% compute scattering transform of first image
x = image_set(:,:,1);
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);
Sx = scat(x, Wop);
Sx_mat = format_scat(Sx);

% store scattering invariants (one row per image)
sx_all = zeros(m, size(Sx_mat, 1));

% compute scattering invariants for each image
for i=1:m
    x = image_set(:,:,i);
    
    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

%% dmaps
W = squareform(pdist(sx_all)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

if corr(V(:,2), L(:,1)) < 0
    V(:,2) = -V(:,2);
end

figure;
plot(L(:,1), V(:,2),'.')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
if print_figures
    print('DMAPS_scat_time_corr',fmt,res)
end

[~, I] = sort(V(:,2));
if print_figures
    figure;
    for i=im_save_idx
        %subplot(subplot_dim1,subplot_dim2,i)

        imshow(uint8(image_set(:,:,I(i))), 'InitialMagnification', 'fit')
        % make green colormap
        cm_green = gray;
        cm_green(:,1) = 0;
        cm_green(:,3) = 0;
        colormap(cm_green)
        axis off
        print(sprintf('dpERK_scat_%d',i),fmt,res)
        clf
    end
end