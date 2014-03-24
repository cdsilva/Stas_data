clear all
%close all

time_data = '../membrane_pictures/large_dataset/time.mat';
profile_data = '../membrane_pictures/large_dataset/dpERK_DL_140313.mat';

res = '-r300';
fmt = '-djpeg';
fontsize = 30;

im_save_dir = 'paper_figures';

%% load membrane lengths

load(time_data);
mem_lengths = length;
clear length;

% remove some younger embryos
ind = setdiff(1:90, [3, 20, 25, 32, 53, 59, 61, 66, 82]);

m = length(ind);

mem_lengths = mem_lengths(ind);

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% load profiles

load(profile_data);

data_red = dpERK_big_feb7;
data_green = Dl_big_feb7;

[m, n] = size(data_red);

nchannels = 3;
dpERK = zeros(m, n, nchannels);
dpERK(:, :, 1) = data_red;
dpERK(:, :, 2) = data_green;

cm_red = gray;
cm_red(:,2) = 0;
cm_red(:,3) = 0;

for i=1:2
    dpERK(:,:,i) = dpERK(:,:,i) / max(max(dpERK(:,:,i))) * 256;
end

% scramble data alignments
dpERK_unaligned = zeros(size(dpERK));
rng(12345);
rand_offsets = zeros(m,1);
dpERK_unaligned(1,:,:) = dpERK(1,:,:);
for i=2:m
    rand_offsets(i) = randi(n);
    dpERK_unaligned(i,:,:) = circshift(dpERK(i,:,:),[0 rand_offsets(i) 0]);
end

%%
figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(uint8(dpERK_unaligned))
xlabel('position (unregistered)', 'fontsize', fontsize)
ylabel('data sample (unordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca,'position',[0.1 0.1 0.9 0.9],'units','normalized')

print(sprintf('%s/unregistered_unordered_1d', im_save_dir), fmt, res);

%% 

idx = 16;

npoints = 500;
[X, Y] = meshgrid(1:npoints, 1:npoints);

X = X - mean(X(:));
Y = Y - mean(Y(:));

R = sqrt(X.^2 + Y.^2);
theta = atan2(Y, X);
idx2 = find(theta < 0);
theta(idx2) = theta(idx2) + 2*pi;

image_test = zeros(npoints, npoints, nchannels);
for j=1:nchannels
    image_test_tmp = zeros(npoints);
    for i=1:n
        idx2 = find(theta > (i-1)/n*2*pi & theta < i/n * 2 * pi & R > 0.7*(npoints/2) & R < 0.8*(npoints/2));
        image_test_tmp(idx2) = dpERK_unaligned(idx, i, j);
    end
    image_test(:,:,j) = image_test_tmp;
end

figure;
set(gcf, 'paperposition',[0 0 9 10])
subplot(10, 1, 1:9)
imshow(uint8(image_test))
axis off
set(gca, 'ydir','normal')

subplot(10, 1, 10)
imshow(uint8(dpERK_unaligned(idx,:,:)));
axis off

print(sprintf('%s/illustrate_1d', im_save_dir), fmt, res);

%%

dim = 2;

[R, W] = compute_pairwise_alignments_1d_color(dpERK_unaligned);

eps = median(W(:));
neigs = 5;

%[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
R_opt = ang_synch(R, dim);

%%

dpERK_aligned = zeros(size(dpERK_unaligned));

for i=1:m
    R_tmp = R_opt(dim*(i-1)+1:dim*i,:);
    theta = atan2(R_tmp(1,2), R_tmp(1,1));
    dpERK_aligned(i,:,:) = circshift(dpERK_unaligned(i,:,:), [0 round(theta/(2*pi)*n) 0]);
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(uint8(dpERK_aligned))
xlabel('position (registered)', 'fontsize', fontsize)
ylabel('data sample (unordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca,'position',[0.1 0.1 0.9 0.9],'units','normalized')

print(sprintf('%s/registered_unordered_1d', im_save_dir), fmt, res);

W2 = zeros(m);
for i=1:m
    imagei = dpERK_aligned(i,:,:);
    for j=1:i-1    
        imagej = dpERK_aligned(j,:,:);
        W2(i, j) = norm(imagei(:)-imagej(:))^2;
        W2(j,i) = W2(i,j);
    end
end
eps2 = median(W2(:));

[V, D] = dmaps(W2, eps2, neigs);

if corr(V(:,2), mem_lengths) < 0
    V(:,2) = -V(:,2);
end

[~, I] = sort(V(:,2));

figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(uint8(dpERK_aligned(I, :,:)))
xlabel('position (registered)', 'fontsize', fontsize)
ylabel('data sample (ordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca,'position',[0.1 0.1 0.9 0.9],'units','normalized')

print(sprintf('%s/registered_ordered_1d', im_save_dir), fmt, res);

ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_angsynch = compute_ranks(V(:,2));


figure;
plot(ranks_from_membranes, ranks_from_angsynch, '.')
xlabel('rank from membrane lengths')
ylabel('rank from synchronization + dmaps')


%%

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

dpERK_aligned = zeros(size(dpERK_unaligned));

for i=1:m
    R_tmp = R_opt(dim*(i-1)+1:dim*i,:);
    theta = atan2(R_tmp(1,2), R_tmp(1,1));
    dpERK_aligned(i,:,:) = circshift(dpERK_unaligned(i,:,:), [0 round(theta/(2*pi)*n) 0]);
end

figure;
imshow(uint8(dpERK_aligned))

idx = find(embed_idx(1,:) == dim+1 & embed_idx(2,:) == 1);

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

[~, I] = sort(embed_coord(:, idx));
figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(uint8(dpERK_aligned(I, :,:)))
colormap(cm_red)
xlabel('position (registered)', 'fontsize', fontsize)
ylabel('data sample (ordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca,'position',[0.1 0.1 0.9 0.9],'units','normalized')
print(sprintf('%s/registered_ordered_vdm_1d', im_save_dir), fmt, res);

ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_vdm = compute_ranks(embed_coord(:,idx));

figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('rank from membrane lengths')
ylabel('rank from vdm')






