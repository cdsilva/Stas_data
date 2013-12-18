clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';

%% load data
load ../dpERK_oct16.mat

%% PCA for each data set

[m, n] = size(dpERK);

dpERK2 = 1+max(max(dpERK)) * ones(2*size(dpERK,1),size(dpERK,2));
dpERK2(2:2:end,:) = dpERK;

jet2 = jet;
jet2(end,:) = [1,1,1];

% plot data
figure;
imagesc(dpERK2)
ylabel('data point (scrambled)')
xlabel('position')
colormap(jet2)
print('data_spaced_scrambled',fmt,res)

% PCA

[V_PCA, D_PCA] = PCA(dpERK, 5);

% mean - center
mean_data = mean(dpERK, 1);
data_PCA = dpERK - ones(m, 1) * mean_data;

if sum(V_PCA(:,1)) < 0
    V_PCA(:,1) = -V_PCA(:,1);
end
coeff = data_PCA * V_PCA(:,1);
[~, I] = sort(coeff);

% data ordered by PCA
dpERK2(2:2:end,:) = dpERK(I,:);
figure;
imagesc(dpERK2)
ylabel('time')
xlabel('position')
colormap(jet2)
print('data_spaced_unscrambled',fmt,res)

coeff = data_PCA * V_PCA;

% plot first two PCA coefficients
figure;
scatter(coeff(:,1),coeff(:,2),200,'b','.')
xlabel('a_{i,1}')
ylabel('a_{i,2}')
axis equal
print('coeff_12_new', fmt, res)

% dmaps
W = squareform(pdist(dpERK));
[V_dmaps, D_dmaps] = dmaps(W, 100, 10);
[~,I] = sort(V_dmaps(:,2));

%p = polyfit(coeff(:,1),coeff(:,2),6);

% plot first two PCA coefficients colored by dmaps
figure;
scatter(coeff(:,1),coeff(:,2),200,V_dmaps(:,2),'.')
hold on
%plot(sort(coeff(:,1)),polyval(p,sort(coeff(:,1))))
coeff = coeff(I,:);
plot(coeff(:,1), smooth(coeff(:,2), 'rlowess'))
xlabel('a_{i,1}')
ylabel('a_{i,2}')
axis equal
print('coeff_12_colored_new', fmt, res)


