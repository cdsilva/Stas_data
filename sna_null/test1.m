clear all 
close all

load Snanull_mar14;

dpERK = [m_dpERK; w_dpERK];
Dl = [m_Dl; w_Dl];
is_mutant = [true(size(m_dpERK, 1), 1); false(size(w_dpERK, 1), 1)];

%%
W = squareform(pdist(m_dpERK));
eps = median(W(:));

[m_V, m_D] = dmaps(W, eps, 10);
[~, m_I] = sort(m_V(:,2));

if sum(m_dpERK(m_I(1))) > sum(m_dpERK(m_I(end)))
    m_V(:,2) = -m_V(:,2);
    [~, m_I] = sort(m_V(:,2));
end

%%
W = squareform(pdist(w_dpERK));
eps = median(W(:));

[w_V, w_D] = dmaps(W, eps, 10);
[~, w_I] = sort(w_V(:,2));

if sum(w_dpERK(w_I(1))) > sum(w_dpERK(w_I(end)))
    w_V(:,2) = -w_V(:,2);
    [~, w_I] = sort(w_V(:,2));
end

%%
figure;
imagesc(m_dpERK(m_I, :));
xlabel('position')
ylabel('ordered by \phi_2')
title('sna null data')

figure;
imagesc(w_dpERK(w_I, :));
xlabel('position')
ylabel('ordered by \phi_2')
title('wt data')

%%
W = squareform(pdist(dpERK));
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

[~, I] = sort(V(:,2));

figure; 
plot(V(is_mutant,2),V(is_mutant,3),'.b')
hold on
plot(V(~is_mutant,2),V(~is_mutant,3),'.r')

m_ind(m_I) = 1:length(m_I);
w_ind(w_I) = 1:length(w_I);

figure; 
scatter(V(is_mutant,2),V(is_mutant,3), 50, m_ind, 's')
hold on
scatter(V(~is_mutant,2),V(~is_mutant,3), 200, w_ind, '.')
legend('sna null','wt','location','best')
xlabel('\phi_2')
ylabel('\phi_3')
colorbar

figure;
imagesc(dpERK(I, :))
