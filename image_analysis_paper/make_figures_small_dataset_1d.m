clear all
close all

time_data = '../membrane_lengths/oct16.mat';

%% load membrane lengths

load(time_data);
mem_lengths = L(:,1);
clear length;

m = length(mem_lengths);
ind = 1:m;

% compute ranks from membrane lengths
ranks_from_membranes = compute_ranks(mem_lengths);

%% load profiles

dpERK = dpERK_raw;

dpERK = dpERK(ind, :);

[m, n] = size(dpERK);

% scramble data alignments
dpERK_unaligned = zeros(size(dpERK));
rng(12345);
rand_offsets = zeros(m,1);
for i=1:m
    rand_offsets(i) = randi(n);
    dpERK_unaligned(i,:) = circshift(dpERK(i,:),[0 rand_offsets(i)]);
end

%%

figure;
imagesc(dpERK_unaligned)
colormap(hot)

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
imagesc(dpERK_aligned)
colormap(hot)

W2 = squareform(pdist(dpERK_aligned));
eps2 = median(W2(:));

[V, D] = dmaps(W2, eps2, neigs);

if corr(V(:,2), mem_lengths) < 0
    V(:,2) = -V(:,2);
end

[~, I] = sort(V(:,2));

figure;
imagesc(dpERK_aligned(I, :))
colormap(hot)

%%

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

dpERK_aligned = zeros(size(dpERK_unaligned));

for i=1:m
    R_tmp = R_opt(dim*(i-1)+1:dim*i,:);
    theta = atan2(R_tmp(1,2), R_tmp(1,1));
    dpERK_aligned(i,:) = circshift(dpERK_unaligned(i,:), [0 round(theta/(2*pi)*n)]);
end

figure;
imagesc(dpERK_aligned)
colormap(hot)

idx = find(embed_idx(1,:) == dim+1 & embed_idx(2,:) == 1);

if corr(mem_lengths, embed_coord(:, idx)) < 0
    embed_coord(:, idx) = -embed_coord(:, idx);
end

[~, I] = sort(embed_coord(:, idx));
figure;
imagesc(dpERK_aligned(I, :))
colormap(hot)

ranks_from_membranes = compute_ranks(mem_lengths);
ranks_from_vdm = compute_ranks(embed_coord(:,idx));


figure;
plot(ranks_from_membranes, ranks_from_vdm, '.')
xlabel('rank from membrane lengths')
ylabel('rank from vdm')






