clear all
close all

%% directories where things are stored
time_data = '../membrane_lengths/oct16.mat';
image_dir = '../membrane_pictures/membrane2/dpERK_staining';

im_save_dir = 'paper_figures2';


%% load membrane lengths

load(time_data);
mem_lengths = L(:,1);

% remove some embryos
ind = setdiff(1:52, [9, 15, 22]);

m = length(ind);

mem_lengths = mem_lengths(ind);

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% image parameters

npixels = 101;

%set image plotting parameters
subplot_dim1 = ceil(sqrt(m));
subplot_dim2 = ceil(m / subplot_dim1);
[Y, X] = meshgrid((subplot_dim1:-1:1)/subplot_dim1, (1:subplot_dim2)/subplot_dim2);

data = zeros(npixels, npixels, m);
data2 = zeros(npixels, npixels, m);
image_set = zeros(npixels, npixels, 3, m, 'uint8');


%% load images

for i=1:m
    % read image
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, ind(i)));
    
    
    % resize image
    im1 = imresize(im1, [npixels npixels]);
    
    
    %store image
    data(:, :, i) = im1(:, :, 2);
    data2(:, :, i) = im1(:, :, 3);
    
    im1 = circshift(im1, [0 0 -1]);
    im1(:, :, 1) = im1(:, :, 1) + im1(:, :, 3);
    im1(:, :, 2) = im1(:, :, 2) + im1(:, :, 3);
    
    image_set(:, :, :, i) = im1;
    
end

%%

r_max = floor(npixels/2);
K = m;

[ I, Eigim, Freqs, Rad_Freqs, W ] = denoise_images(data, r_max,  K);

%%

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [16 16])
set(gcf, 'paperposition',[0 0 16 16 ]);
for i=1:m
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(uint8(data(:, :, i)));
end
% saveas(gcf,sprintf('%s/FBsPCA_images', im_save_dir), 'pdf')


figure;
for i=1:m
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(uint8(I(:, :, i)));
end


%%
npixels=size(data, 1);
m=size(data, 3);
N=floor(npixels/2);
trunc = 200; % maximum number of FBsPCA leading componnents chosen.

[x, y]=meshgrid(-N:N, -N:N);    %generate square grids
r=sqrt(x.^2+y.^2);      %compute the radius at each grid point
test=reshape(data, npixels^2, m);
test=test(r>r_max, :);
noise_variance=var(test(:)); %Estimate the noise variance from ear region.

%Use Fourier Bessel steerable PCA to denoise and compress the original data
%set.
[ U, D, freqs, rad_freqs, Mean ] = FB_SPCA(data, r_max);
[ UU, Freqs, ~, W ] = FBSPCA_MP_rankEst(m, U, D, freqs, rad_freqs, max(noise_variance, D(trunc)));
%Cutoff is determined by the maximum of the noise_variance and the 200th
%eigenvalues.
fprintf('\nThe number of components is %d \n', size(Freqs, 1));

% Using Linear Wiener Filter to denoise:
[ Coeff ] = WF_FBSPCA( data, Mean, r_max, [0, 0], UU, Freqs, W, 1);

% Coeff(Freqs==0, :)=Coeff(Freqs==0, :)/sqrt(2);
% for i=1:m
%     Coeff(:, i)=Coeff(:, i)/norm(Coeff(:, i));
% end;
% Coeff(Freqs==0, :)=Coeff(Freqs==0, :)*sqrt(2);

[ Coeff_b, Coeff_b_r, toc_bispec ] = Bispec_2Drot( Coeff, Freqs );


%%
npixels=size(data2, 1);
m=size(data2, 3);
N=floor(npixels/2);
trunc = 200; % maximum number of FBsPCA leading componnents chosen.

[x, y]=meshgrid(-N:N, -N:N);    %generate square grids
r=sqrt(x.^2+y.^2);      %compute the radius at each grid point
test=reshape(data2, npixels^2, m);
test=test(r>r_max, :);
noise_variance=var(test(:)); %Estimate the noise variance from ear region.

%Use Fourier Bessel steerable PCA to denoise and compress the original data
%set.
[ U, D, freqs, rad_freqs, Mean ] = FB_SPCA(data2, r_max);
[ UU, Freqs, ~, W ] = FBSPCA_MP_rankEst(m, U, D, freqs, rad_freqs, max(noise_variance, D(trunc)));
%Cutoff is determined by the maximum of the noise_variance and the 200th
%eigenvalues.
fprintf('\nThe number of components is %d \n', size(Freqs, 1));

% Using Linear Wiener Filter to denoise:
[ Coeff2 ] = WF_FBSPCA( data2, Mean, r_max, [0, 0], UU, Freqs, W, 1);

% Coeff(Freqs==0, :)=Coeff(Freqs==0, :)/sqrt(2);
% for i=1:m
%     Coeff(:, i)=Coeff(:, i)/norm(Coeff(:, i));
% end;
% Coeff(Freqs==0, :)=Coeff(Freqs==0, :)*sqrt(2);

[ Coeff_b2, Coeff_b_r2, toc_bispec ] = Bispec_2Drot( Coeff2, Freqs );

%%
% W_dmaps = zeros(m);
% for i=1:m
%     for j=1:i-1
%         W_dmaps(i, j) = norm(Coeff_b(3,i) - Coeff_b(3,j)).^2;
%     end
% end
% W_dmaps = W_dmaps + W_dmaps';
W_dmaps = squareform(pdist([real(Coeff_b') imag(Coeff_b') real(Coeff_b2') imag(Coeff_b2')])).^2;
eps_dmaps = median(W_dmaps(:))/4;

[V, D] = dmaps(W_dmaps, eps_dmaps, m);

figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [8 8])
set(gcf, 'paperposition',[0 0 8 8]);
plot(diag(abs(D)),'.')
xlabel('k')
ylabel('|\lambda_k|')
saveas(gcf,sprintf('%s/FBsPCA_evals', im_save_dir), 'pdf')

figure;
imshow(W_dmaps(ranks_from_membranes,ranks_from_membranes))

ncomp = 6;
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [24 24/ncomp*2])
set(gcf, 'paperposition',[0 0 24 24/ncomp*2]);
for i=2:ncomp+1
    subplot(2, ncomp, i-1)
    plot(V(:,i),mem_lengths,'.')
    xlabel(sprintf('\\phi_%d',i));
    ylabel('mem. length')
    axis([-1.1 1.1 -inf inf])
    
    subplot(2, ncomp, ncomp+(i-1))
    [~, idx] = max(abs(V(:,i)));
    imshow(image_set(:, :, :, idx), 'border', 'tight')
    
end
saveas(gcf,sprintf('%s/FBsPCA_dmaps', im_save_dir), 'pdf')


%%
% W_dmaps = zeros(m);
% for i=1:m
%     for j=1:i-1
%         W_dmaps(i, j) = sum((abs(Coeff(:,i)) - abs(Coeff(:,j))).^2);
%     end
% end
% W_dmaps = W_dmaps + W_dmaps';
% eps_dmaps = median(W_dmaps(:));
%
% [V, D] = dmaps(W_dmaps, eps_dmaps, 10);
%
% figure;
% plot(V(:,2),mem_lengths,'.')

[~, I] = sort(V(:,2));
figure;
for i=1:m
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(uint8(data(:, :, I(i))));
end

%%
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [16 16])
set(gcf, 'paperposition',[0 0 16 16 ]);
for i=1:m
    make_subplot(subplot_dim1, subplot_dim2, 0.01, i);
    imshow(uint8(data(:, :, i)));
end



