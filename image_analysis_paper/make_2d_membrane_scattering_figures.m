
% pictures are in image channel 2
image_channel = 2;
% number of pixels (subsample images)
npixels = 100;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);

% store images in 3d array
image_set = zeros(npixels, npixels, m);

% load in images
figure;
for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', dpERK_membrane_dir, i));

    % resize image
    im1 = imresize(im1, [npixels npixels]);

    % extract relevent color from image
    im1 = im1(:,:,image_channel);

    % plot image
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(im1)

    %store image
    im1 = double(im1);
    image_set(:,:,i) = im1;
end

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
ylabel(sprintf('\\phi_%d',idx))
if print_figures
    print('DMAPS_membrane_scat_time_corr',fmt,res)
end

[~, I] = sort(V(:,idx));
if print_figures
    figure;
    for i=1:m
        %subplot(subplot_dim1,subplot_dim2,i)

        imshow(uint8(image_set(:,:,I(i))), 'InitialMagnification', 'fit')
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