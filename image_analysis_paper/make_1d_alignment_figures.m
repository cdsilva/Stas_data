addpath ../membrane_pictures/synchron_SO2

%% plot scrambled data
figure;
imagesc(dpERK_unaligned)
xlabel('position (unaligned)')
ylabel('data point (unordered)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_unaligned_unordered',fmt, res)
end

%% calculate pairwise alignments
tic
[H, W, theta] = align_data(dpERK_unaligned);

%% angular synchronization
R_opt = ang_synch(H, 2);

toc
dpERK_aligned = zeros(size(dpERK_unaligned));
R_est = R_opt(1:2, 1:2);
theta_est = atan2(R_est(2,1),R_est(1,1));
shift_est0 = round(theta_est/(2*pi)*n) + rand_offsets(1);
for i=1:m
    R_est = R_opt(2*i-1:2*i, 1:2);
    theta_est = atan2(R_est(2,1),R_est(1,1));
    shift_est = round(theta_est/(2*pi)*n) - shift_est0;
    dpERK_aligned(i,:) = circshift(dpERK_unaligned(i,:),[0 shift_est]);
end

%% plot alignment with angular synchronization
figure;
imagesc(dpERK_aligned)
xlabel('position (aligned using angular synchronization)')
ylabel('data point (unordered)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_aligned_unordered',fmt, res)
end

%% dmaps following angular synchronization
W = squareform(pdist(dpERK_aligned)).^2;
eps_dmaps = median(W(:));
[V_dmaps, D_dmaps] = dmaps(W, eps_dmaps, 10);

if corr(V_dmaps(:,2), L(:,1)) < 0
    V_dmaps(:,2) = -V_dmaps(:,2);
end

%% plot dmaps ordering
[~, I] = sort(V_dmaps(:,2));
figure;
imagesc(dpERK_aligned(I,:))
ylabel('data point (ordered using DMAPS)')
xlabel('position (aligned using angular synchronization)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_ordered_angsynch',fmt, res)
end

%% plot correlation
figure;
plot(L(:,1), V_dmaps(:,2), '.')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
if print_figures
    print('angsynch_time_corr',fmt, res)
end

fprintf('Angular synchronization Spearman coeff: %2.4f \n', corr(L(:,1),V_dmaps(:,2), 'type','spearman'));

    
%% VDM
% calculate pairwise alignments
tic
[H, W, theta] = align_data(dpERK_unaligned);

eps = median(W(:));
neigs = 6;

[R_opt, embed_coord, embed_idx, D] = vdm(H, W, eps, neigs);
toc

data2 = zeros(size(dpERK_unaligned));
R_est = R_opt(1:2,1:2);
theta_est = atan2(R_est(2,1),R_est(1,1));
shift_est0 = round(theta_est/(2*pi)*n) + rand_offsets(1);
for i=1:m
    R_est = R_opt(2*i-1:2*i,:);
    theta_est = atan2(R_est(2,1),R_est(1,1));
    shift_est = round(theta_est/(2*pi)*n) - shift_est0;
    data2(i,:) = circshift(dpERK_unaligned(i,:),[0 shift_est]);
end

%% plot aligned data
figure;
imagesc(data2)
xlabel('position (aligned using VDM)')
ylabel('scrambled data')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_aligned_vdm',fmt, res)
end

% figure;
% for i=1:5
%     subplot(3,2,i)
%     plot(L(:,1), embed_coord(:,i), '.')
% end

coord_idx = 4;
if corr(embed_coord(:,coord_idx), L(:,1)) < 0
    embed_coord(:,coord_idx) = -embed_coord(:,coord_idx);
end
[~, idx] = sort(embed_coord(:,coord_idx));

%% plot aligned and ordered data
figure;
imagesc(data2(idx,:))
ylabel('data point (ordered using VDM)')
xlabel('position (aligned using VDM)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_ordered_vdm',fmt, res)
end

%% plot correlation
figure;
plot(L(:,1), embed_coord(:, coord_idx), '.')
xlabel('membrane thickness')
ylabel(sprintf('$\\langle \\phi_%d, \\phi_%d \\rangle$', embed_idx(1, coord_idx), embed_idx(2, coord_idx)),'interpreter','latex')
if print_figures
    print('VDM_time_corr',fmt, res)
end

fprintf('VDM (1-D) Spearman coeff: %2.4f \n', corr(L(:,1),embed_coord(:, coord_idx), 'type','spearman'));

rmpath ../membrane_pictures/synchron_SO2