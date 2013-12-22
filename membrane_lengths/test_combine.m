clear all
close all

load curve_data.mat
x = interp1(xaxis,time,x);
y = interp1(yaxis,furrow_len,y);
xmin = interp1(xaxis,time,xmin);
ymin = interp1(yaxis,furrow_len,ymin);
xmax = interp1(xaxis,time,xmax);
ymax = interp1(yaxis,furrow_len,ymax);

files = {'oct16.mat'; 'feb11.mat'; 'mar15.mat'};
nfiles = size(files,1);

new_cm = jet(10);
new_cm = new_cm(2:end-1,:);

%% combine all data; subtract out min average
data = [];
times = [];
data_ind = [];
for i=1:nfiles
    load(files{i});
    min(mean(dpERK_raw,2))
    mean(mean(dpERK_raw,2))
    data = [data; dpERK_raw-mean(mean(dpERK_raw,2))];
    % data = [data; dpERK_raw];
    times =  [times; interp1(y, x, L(:,1), [], 'extrap')];
    data_ind = [data_ind; i*ones(size(L(:,1)))];
end

[~, I_true] = sort(times);

%% order by PCA
[I, coeff, V_PCA, D_PCA] = unscramble_pca(data);

figure;
imagesc(data(I,:))

figure;
imagesc(data(I_true,:))

figure; 
scatter(1:length(I),times(I),50, data_ind(I),'.')
colormap(new_cm)

figure;
scatter(times, -coeff(:,1),50, data_ind(I),'.')
colormap(new_cm)
xlabel('time (from membrane length)')
ylabel('PCA proj coeff')

%% regression
ind = find(data_ind == 1);
% X = [ones(length(data_ind),1) coeff(:, 1:2) coeff(:,1:2).^2];
X = [ones(length(data_ind),1) coeff(:, 1)];
B = regress(times(ind), X(ind,:));

ind2 = find(data_ind == 2);
times_est = X(ind2,:) * B;

figure;
plot(times(ind2), times_est, '.')
hold on
plot([0 65], [0 65], 'r')

%% order by DMAPS
% eps = logspace(3, 7, 10);
% sumA = zeros(size(eps));
% W = squareform(pdist(data)).^2;
% for i=1:length(eps)
%     sumA(i) = sum(sum(exp(-W/eps(i))));
% end
% 
% figure; 
% loglog(eps, sumA)
% 
% eps = 5e4;
% [I, V_DMAPS, D_DMAPS] = unscramble_dmaps(data, eps);
% 
% figure;
% imagesc(data(I,:))
% 
% figure; 
% scatter(1:length(I),times(I),50, data_ind(I),'.')

