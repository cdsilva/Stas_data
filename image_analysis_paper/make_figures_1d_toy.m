clear all
close all


res = '-r300';
fmt = '-djpeg';
fontsize = 20;

im_save_dir = 'paper_figures';

rng(12345);

%% load membrane lengths

m = 100;

mem_lengths = randperm(m)';

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% load profiles

n = 100;

dpERK = zeros(m, n);

sigma = 100;
for i=1:m
    dpERK(i, :) = mem_lengths(i)*exp(-((1:n)-(n/2)).^2/sigma);
end
sigma_noise = 5;
dpERK = dpERK + sigma_noise*randn(size(dpERK));

cm_red = gray;
cm_red(:,2) = 0;
cm_red(:,3) = 0;

dpERK = dpERK / max(dpERK(:)) * size(cm_red, 1);

% scramble data alignments
dpERK_unaligned = zeros(size(dpERK));
rand_offsets = zeros(m,1);
dpERK_unaligned(1,:) = dpERK(1,:);
for i=2:m
    rand_offsets(i) = randi(n);
    dpERK_unaligned(i,:) = circshift(dpERK(i,:),[0 rand_offsets(i)]);
end


%%
figure;
set(gcf, 'paperposition',[0 0 8 8])
image(dpERK_unaligned)
colormap(cm_red)
xlabel('position (unregistered)', 'fontsize', fontsize)
ylabel('data sample (unordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print(sprintf('%s/unregistered_unordered_toy', im_save_dir), fmt, res);

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

image_test = zeros(npoints);
for i=1:n
    idx2 = find(theta > (i-1)/n*2*pi & theta < i/n * 2 * pi & R > 0.7*(npoints/2) & R < 0.8*(npoints/2));
    image_test(idx2) = dpERK_unaligned(idx, i);
end

figure;
set(gcf, 'paperposition',[0 0 9 10])
subplot(10, 1, 1:9)
image(image_test)
axis off
set(gca, 'ydir','normal')
colormap(cm_red)

subplot(10, 1, 10)
image(dpERK_unaligned(idx,:));
axis off
colormap(cm_red)

print(sprintf('%s/illustrate_toy', im_save_dir), fmt, res);

%%

dim = 2;

[R, W] = compute_pairwise_alignments_1d(dpERK_unaligned);

eps = median(W(:));
neigs = 5;

%[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
R_opt = ang_synch(R, dim);

%%

dpERK_aligned = zeros(size(dpERK_unaligned));

for i=1:m
    R_tmp = R_opt(dim*(i-1)+1:dim*i,:);
    theta = atan2(R_tmp(1,2), R_tmp(1,1));
    dpERK_aligned(i,:) = circshift(dpERK_unaligned(i,:), [0 round(theta/(2*pi)*n)]);
end

figure;
set(gcf, 'paperposition',[0 0 8 8])
image(dpERK_aligned)
colormap(cm_red)
xlabel('position (registered)', 'fontsize', fontsize)
ylabel('data sample (unordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print(sprintf('%s/registered_unordered_toy', im_save_dir), fmt, res);

W2 = squareform(pdist(dpERK_aligned));
eps2 = median(W2(:));

[V, D] = dmaps(W2, eps2, neigs);

if corr(V(:,2), mem_lengths) < 0
    V(:,2) = -V(:,2);
end

[~, I] = sort(V(:,2));

figure;
set(gcf, 'paperposition',[0 0 8 8])
image(dpERK_aligned(I, :))
colormap(cm_red)
xlabel('position (registered)', 'fontsize', fontsize)
ylabel('data sample (ordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print(sprintf('%s/registered_ordered_toy', im_save_dir), fmt, res);

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
    dpERK_aligned(i,:) = circshift(dpERK_unaligned(i,:), [0 round(theta/(2*pi)*n)]);
end

figure;
image(dpERK_aligned)
colormap(cm_red)

idx = find(embed_idx(1,:) == dim+1 & embed_idx(2,:) == 1);

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

[~, I] = sort(embed_coord(:, idx));
figure;
set(gcf, 'paperposition',[0 0 8 8])
image(dpERK_aligned(I, :))
colormap(cm_red)
xlabel('position (registered)', 'fontsize', fontsize)
ylabel('data sample (ordered)', 'fontsize', fontsize)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print(sprintf('%s/registered_ordered_vdm_toy', im_save_dir), fmt, res);

ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_vdm = compute_ranks(embed_coord(:,idx));


figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('rank from membrane lengths')
ylabel('rank from vdm')






