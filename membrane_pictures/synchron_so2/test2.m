clear all
close all

res = '-r300';
fmt = '-djpeg';

%% load data
load('../../dpERK_oct16.mat');
load('../../membrane_lengths/oct16.mat')
mem_lengths = L(:,1);

data0 = dpERK;
[m, n] = size(data0);

%% dmaps

W = squareform(pdist(data0)).^2;
eps = median(median(W));

[V, D] = dmaps(W, eps, 10);

[~, ind] = sort(V(:,2));

figure;
imagesc(data0)
xlabel('position')
ylabel('data index (unordered)')
%print('prealigned_unordered',fmt, res)

figure;
imagesc(data0(ind,:))
ylabel('data index (ordered)')
%print('prealigned_ordered',fmt, res)

figure;
if sign(mem_lengths(2)-mem_lengths(1)) ~= sign(V(2,2)-V(1,2))
    V(:,2) = -V(:,2);
end
plot(mem_lengths, V(:,2),'.')
xlabel('membrane length')
ylabel('\phi_2')
%print('dmaps_membrane_corr',fmt, res)

%% scramble data alignments
rng(12345);

rand_offsets = zeros(m,1);
for i=1:m
    rand_offsets(i) = randi(n);
    data0(i,:) = circshift(data0(i,:),[0 rand_offsets(i)]);
end

figure;
imagesc(data0)
xlabel('position (unaligned)')
ylabel('data index (unordered)')
%print('unaligned_unordered',fmt, res)

%% calculate pairwise alignments

[R, W, theta] = align_data(data0);

%% angular synchronization
[V, D] = eigs(R, 2);

data = zeros(size(data0));
[u, s, v] = svd(V(1:2,1:2));
R_est = u*v';
theta_est = atan2(R_est(2,1),R_est(1,1));
shift_est0 = round(theta_est/(2*pi)*n) + rand_offsets(1);
for i=1:m
    [u, s, v] = svd(V(2*i-1:2*i,:));
    R_est = u*v';
    theta_est = atan2(R_est(2,1),R_est(1,1));
    shift_est = round(theta_est/(2*pi)*n) - shift_est0;
    data(i,:) = circshift(data0(i,:),[0 shift_est]);
end

figure;
imagesc(data)
xlabel('position')
ylabel('data index (unordered)')
%print('aligned_unordered',fmt, res)

W_dmaps = squareform(pdist(data)).^2;
eps_dmaps = median(median(W_dmaps));
[V_dmaps, D_dmaps] = dmaps(W_dmaps, eps_dmaps, 10);
[~, idx] = sort(V_dmaps(:,2));
figure;
imagesc(data(idx,:))
xlabel('position')
ylabel('data index (ordered)')
%print('aligned_ordered',fmt, res)
    
%% vector dmaps

eps = median(W(:));
neigs = 6;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);

data2 = zeros(size(data0));
R_est = R_opt(1:2,1:2);
theta_est = atan2(R_est(2,1),R_est(1,1));
shift_est0 = round(theta_est/(2*pi)*n) + rand_offsets(1);
for i=1:m
    R_est = R_opt(2*i-1:2*i,:);
    theta_est = atan2(R_est(2,1),R_est(1,1));
    shift_est = round(theta_est/(2*pi)*n) - shift_est0;
    data2(i,:) = circshift(data0(i,:),[0 shift_est]);
end

figure;
imagesc(data2)
xlabel('position (unaligned)')
ylabel('scrambled data')

coord_idx = 2;
[~, idx] = sort(embed_coord(:,coord_idx));
figure;
imagesc(data2(idx,:))


    
    