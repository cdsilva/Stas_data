close all

%% do computations
tic

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
    %[~, thres] = edge(uint8(image_set_membrane(:,:,i)), 'log');
    %x = edge(uint8(image_set_membrane(:,:,i)), 'log');
    x = image_set_membrane(:,:,i);
    Sx = scat(x, Wop);

    sx_all(i,:) = mean(mean(format_scat(Sx),2),3)';
end

% dmaps

W = squareform(pdist(sx_all)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 20);

runtime = toc

% for i=3:20;
%     figure;
%     plot(V(:,2),V(:,i),'.')
% end
% 
% for i=2:20
%     corr(V(:,i),L(:,1));
% end

%idx = 8;
[~, idx] = max(abs(corr(V, L(:,1))));

if corr(V(:,idx), L(:,1)) < 0
    V(:,idx) = -V(:,idx);
end

[~, I] = sort(V(:,idx));

%% make plots
figure;
plot(L(:,1), V(:,idx),'.')
hold on
plot(L(I(im_save_idx),1), V(I(im_save_idx), idx), 'o')
xlabel('membrane thickness')
ylabel(sprintf('$\\phi_%d$',idx),'interpreter','latex')
if print_figures
    print('DMAPS_membrane_scat_time_corr',fmt,res)
end

fprintf('Scattering transform (membrane) Spearman coeff: %2.4f \n', corr(L(:,1), V(:,idx), 'type','spearman'));

if print_figures
    figure;
    for i=im_save_idx
        set(gcf, 'paperposition',[0 0 8 8])
        imshow(logical(image_set_membrane(:,:,I(i))), 'InitialMagnification', 'fit')
        set(gca,'position',[0 0 1 1],'units','normalized')

        axis off
        print(sprintf('membrane_scat_%d',i),fmt,res)
        clf
        
        imshow(uint8(image_set_membrane_raw(:,:,I(i))), 'InitialMagnification', 'fit')
        set(gca,'position',[0 0 1 1],'units','normalized')

        %imshow(imadjust(uint8(image_set(:,:,I(i)))), 'InitialMagnification', 'fit')
        % make green colormap
        cm_green = gray;
        cm_green(:,1) = 0;
        cm_green(:,3) = 0;
        colormap(cm_green)
        axis off
        print(sprintf('membrane_scat_raw_%d',i),fmt,res)
        clf
    end
end

save('2d_membrane_scattering_figures.mat');
