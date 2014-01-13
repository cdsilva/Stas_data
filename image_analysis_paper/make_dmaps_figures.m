%% distance matrix
tic
W = squareform(pdist(dpERK));
eps = median(W(:));
[V_dmaps, D_dmaps] = dmaps(W, eps, 10);
toc

%% flip sign if necessary
if corr(V_dmaps(:,2), L(:,1)) < 0
    V_dmaps(:,2) = -V_dmaps(:,2);
end

%% order and plot
[~,I] = sort(V_dmaps(:,2));

figure;
imagesc(dpERK(I,:))
ylabel('data point (ordered using DMAPS)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_ordered_DMAPS',fmt,res)
end

%% correlation
figure;
plot(L(:,1),V_dmaps(:,2),'.')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
if print_figures
    print('DMAPS_time_corr',fmt,res)
end

fprintf('DMAPS Spearman coeff: %2.4f \n', corr(L(:,1), V_dmaps(:,2), 'type','spearman'));

save('dmaps_figures.mat');
