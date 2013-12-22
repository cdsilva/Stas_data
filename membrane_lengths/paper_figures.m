clear all
close all

%% define files
files = {'oct16.mat'; 'feb11.mat'; 'mar15.mat'};
nfiles = size(files,1);
color = 'brg';
nboot = 1000;

%% load data from calibration curve
load curve_data.mat
x = interp1(xaxis,time,x);
y = interp1(yaxis,furrow_len,y);
xmin = interp1(xaxis,time,xmin);
ymin = interp1(yaxis,furrow_len,ymin);
xmax = interp1(xaxis,time,xmax);
ymax = interp1(yaxis,furrow_len,ymax);

%% do analysis with first file
load(files{1});

%% plot raw data
figure;
imagesc(dpERK_raw)
xlabel('position')
ylabel('unordered')

%% plot data ordered by membrane length
figure;
imagesc(dpERK_raw(L(:,2),:))
xlabel('position')
ylabel('order by membrane length')

%% denoise data using PCA
[I, coeff, V_PCA, D_PCA] = unscramble_pca(dpERK_raw, true);
data_mean = mean(dpERK_sort);
dpERK_denoise = coeff(:,1) * V_PCA(:,1)' + ones(size(coeff,1),1) * data_mean;

figure;
imagesc(dpERK_denoise(L(:,2),:))
xlabel('position')
ylabel('ordered by membrane length; denoised using PCA')

%% PCA schematic
num_cells = 20;

figure;
imagesc(dpERK_raw(1:num_cells,:))
%axis equal
axis off

figure;
imagesc(coeff(1:num_cells,1))
axis equal
axis off

figure;
imagesc(V_PCA(:,1)')
curr_daspect = daspect;
curr_daspect(1) = curr_daspect(1)/num_cells;
daspect(curr_daspect)
axis off

%% plot modes
figure;
plot(V_PCA(:,1))
xlabel('position')
ylabel('first PCA mode')

figure;
if V_PCA(1,2) < 0
    plot(-V_PCA(:,2))
else
    plot(V_PCA(:,2))
end
xlabel('position')
ylabel('second PCA mode')

figure;
plot(V_PCA(:,3))
xlabel('position')
ylabel('third PCA mode')

figure;
plot(diag(D_PCA),'.')
ylabel('PCA eigenvalue')
set(gca,'xlim',[0 10])

%% plot PCA coeff vs. time
% time = interp1(y, x, L(:,1), [], 'extrap');
% p = polyfit(time, coeff(:,1), 3);
%
% figure;
% plot(time,coeff(:,1),'.')
% hold on
% plot(sort(time), polyval(p,sort(time)))
% xlabel('time from membrane length')
% ylabel('first PCA projection coefficient')

%% compare data sets PCA coeff vs. time
% figure;
% h = zeros(3,1);
% for i=1:nfiles
%     load(files{i,1})
%
%     dpERK_raw = dpERK_raw - min(mean(dpERK_raw,2));
%     [I, coeff, V_PCA, D_PCA] = unscramble_pca(dpERK_raw, true);
%
%     time = interp1(y, x, L(:,1), [], 'extrap');
%     p = polyfit(time, coeff(:,1), 3);
%
%     %     power_law_params = fminsearch(@(params) power_law(params, time, coeff(:,1)+mean(dpERK_raw)*V_PCA(:,1)), [0.1; 3]);
%     %
%     %     %plot(time,coeff(:,1)+mean(dpERK_raw)*V_PCA(:,1),'.', 'color', color(i))
%     %     plot(time,dpERK_raw*V_PCA(:,1),'.', 'color', color(i))
%     %     hold on
%     %     plot(sort(time), power_law_params(1)*sort(time).^power_law_params(2),'color', color(i))
%
%     coeff0 = coeff(:,1)+mean(dpERK_raw)*V_PCA(:,1);
%
%     [power_law_ci, power_law_data] = bootci(100,@power_law_params,time, coeff0);
%
%     power_law_means = power_law_params(time, coeff0);
%     power_law_ci
%     %plot(time,coeff(:,1)+mean(dpERK_raw)*V_PCA(:,1),'.', 'color', color(i))
%     h(i) = plot(time,coeff0,'.', 'color', color(i));
%     hold on
%     plot(sort(time), power_law(power_law_means, sort(time)),'color', color(i))
%     %plot(sort(time), polyval(p,sort(time)),'color', color(i))
%
% end
% xlabel('time from membrane length')
% ylabel('first PCA projection coefficient')
% legend(h,{'data1','data2','data3'},'location','best')

%% compare data sets PCA modes


figure;
for i=1:nfiles
    load(files{i,1})
    
    [m, n] = size(dpERK_raw);
    
    [bootstat_ci, bootstat] = bootci(nboot, {@bootfun_PCA, dpERK_raw});
    
    errorbar(i:3:size(bootstat,2),mean(bootstat(:,i:3:end)),mean(bootstat(:,i:3:end)) - bootstat_ci(1,i:3:end), bootstat_ci(2,i:3:end) - mean(bootstat(:,i:3:end)),'color',color(i))
    
    hold on
    
end
xlabel('position')
ylabel('first PCA mode')
legend('data1','data2','data3','location','best')

%% make master mode
total_data = [];
time = [];
data_ind = [];
for i=1:nfiles
    load(files{i,1})
    
    total_data = [total_data; dpERK_raw-min(mean(dpERK_raw,2))];
    time = [time; interp1(y, x, L(:,1), [], 'extrap')];
    data_ind = [data_ind; i*ones(size(dpERK_raw,1),1)];
end
[bootstat_ci, bootstat] = bootci(nboot, {@bootfun_PCA, total_data});

master_V = mean(bootstat);
[I, coeff, V_PCA, D_PCA] = unscramble_pca(total_data, true);

variance_captured = D_PCA(1,1) / sum(diag(D_PCA))

figure;
errorbar(1:size(bootstat,2),master_V,master_V - bootstat_ci(1,:), bootstat_ci(2,:) -master_V)
xlabel('position')
ylabel('first PCA mode')

%% compare data sets PCA coeff vs. time using master mode
coeff0 = total_data * master_V';
[power_law_ci, power_law_data] = bootci(nboot,@power_law_params,time, coeff0);
power_law_means = power_law_params(time, coeff0);
power_law_ci

figure;
h = zeros(3,1);
for i=1:nfiles
    ind = find(data_ind == i);
    
    h(i) = plot(time(ind),total_data(ind,:) * master_V','.', 'color', color(i));
    hold on
end
plot(sort(time), power_law(power_law_means, sort(time)),'k')

xlabel('time from membrane length')
ylabel('first PCA projection coefficient')
legend(h,{'data1','data2','data3'},'location','best')

%% order using master mode
for i=1:nfiles
    figure;
    load(files{i,1})
    
    [~, I] = sort(dpERK_raw * master_V');
    imagesc(dpERK_raw(I,:))
    xlabel('position')
    ylabel('ordered using PCA')
    set(gca, 'ydir','normal')
end