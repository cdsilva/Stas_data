clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';

%% load data, create arrays to store data
ndata = 3;
data = cell(1,ndata);

load ../dpERK_oct16.mat
data{1,1} = dpERK;

load ../dpERK_jan10.mat
data{1,2} = dpERK_0110;

load ../dpERK_jan17.mat
data{1,3} = dpERK_0117;

c = 'brg';

first_mode = zeros(size(data{1,1},2),ndata);

%% PCA for each data set
for i=1:ndata
    dpERK = data{1,i};
    
    [m, n] = size(dpERK);
    
    % plot data
    figure;
    imagesc(dpERK)
    ylabel('data point (scrambled)')
    xlabel('position')
    print(strcat(['data' num2str(i) '_scrambled']),fmt,res)
    
    % PCA
    
    [V_PCA, D_PCA] = PCA(dpERK, 5);
    
    % mean - center
    mean_data = mean(dpERK, 1);
    data_PCA = dpERK - ones(m, 1) * mean_data;
    
    if sum(V_PCA(:,1)) < 0
        V_PCA(:,1) = -V_PCA(:,1);
    end
    first_mode(:,i) = V_PCA(:,1);
    coeff = data_PCA * V_PCA(:,1);
    [~, I] = sort(coeff);
    
    % plot PCA results
    
    % eigenvalue spectrum
    figure;
    plot(diag(D_PCA),'.')
    xlabel('j')
    ylabel('\lambda_j')
    xlim([0 5])
    print(strcat(['PCA_spectrum' num2str(i)]),fmt,res)
    
    % data ordered by PCA
    figure;
    imagesc(dpERK(I,:))
    ylabel('time')
    xlabel('position')
    print(strcat(['data' num2str(i) '_unscrambled']),fmt,res)
    
    % first PCA mode
    figure;
    plot(V_PCA(:,1),'linewidth',4)
    xlabel('position')
    ylabel('first PCA mode, v_1')
    print(strcat(['PCA_mode' num2str(i)]),fmt,res)
    
    % DMAPS
    W = squareform(pdist(dpERK));
    [V_dmaps, D_dmaps] = dmaps(W, 100, 10);
    [~,I] = sort(V_dmaps(:,2));
    
    % plot data ordered with dmaps
    figure;
    imagesc(dpERK(I,:))
    ylabel('time')
    xlabel('position')
    print(strcat(['data' num2str(i) '_unscrambled_dmaps']),fmt,res)
    
end

%% do PCA and show how we sort coefficients
dpERK = data{1,1};
[m, n] = size(dpERK);

[V_PCA, D_PCA] = PCA(dpERK, 5);

% mean - center
mean_data = mean(dpERK, 1);
data_PCA = dpERK - ones(m, 1) * mean_data;

if sum(V_PCA(:,1)) < 0
    V_PCA(:,1) = -V_PCA(:,1);
end
if sum(V_PCA(:,2)) < 0
    V_PCA(:,2) = -V_PCA(:,2);
end
first_mode(:,i) = V_PCA(:,1);
coeff = data_PCA * V_PCA(:,1:2);
[~,I] = sort(coeff(:,1));

figure;
%plot(coeff(:,1),1:m,'.')
barh(coeff(:,1)-min(coeff(:,1)))
set(gca, 'ydir','reverse')
xlabel('a_{i,1}')
ylabel('i')
ylim([-inf inf])
print('coeff_scrambled', fmt, res)

figure;
plot(coeff(I,1),1:m,'.')
set(gca, 'ydir','reverse')
xlabel('a_{i,1}')
ylabel('i')
print('coeff_unscrambled', fmt, res)

%% motivate dmaps
% plot first two PCA coefficients
figure;
scatter(coeff(:,1),coeff(:,2),200,'b','.')
xlabel('a_{i,1}')
ylabel('a_{i,2}')
axis equal
print('coeff_12', fmt, res)

% dmaps
W = squareform(pdist(dpERK));
[V_dmaps, D_dmaps] = dmaps(W, 100, 10);
[~,I] = sort(V_dmaps(:,2));

% plot first two PCA coefficients colored by dmaps
figure;
scatter(coeff(:,1),coeff(:,2),200,V_dmaps(:,2),'.')
xlabel('a_{i,1}')
ylabel('a_{i,2}')
axis equal
print('coeff_12_colored', fmt, res)
