%% PCA
[V_PCA, D_PCA] = PCA(dpERK, 5);

%% mean - center
mean_data = mean(dpERK, 1);
data_PCA = dpERK - ones(m, 1) * mean_data;

%% flip sign if necessary
if sum(V_PCA(:,1)) < 0
    V_PCA(:,1) = -V_PCA(:,1);
end
coeff = data_PCA * V_PCA;
[~, I] = sort(coeff(:,1));

%% data ordered by PCA
figure;
imagesc(dpERK(I,:))
ylabel('data point (ordered using PCA)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
if print_figures
    print('data_ordered_PCA',fmt,res)
end

%% first PCA mode
figure;
plot(V_PCA(:,1))
xlabel('position')
ylabel('$\psi_1$','interpreter','latex')
if print_figures
    print('PCA_mode',fmt,res)
end

%% correlation
figure;
plot(L(:,1),coeff(:,1),'.')
xlabel('membrane thickness')
ylabel('$\langle x_i, \psi_1 \rangle$','interpreter','latex')
if print_figures
    print('PCA_time_corr',fmt,res)
end

%% plot first two PCA coefficients
figure;
scatter(coeff(:,1),coeff(:,2),200,'b','.')
xlabel('$\langle x_i, \psi_1 \rangle$','interpreter','latex')
ylabel('$\langle x_i, \psi_2 \rangle$','interpreter','latex')
axis equal
if print_figures
    print('coeff_12_new', fmt, res)
end