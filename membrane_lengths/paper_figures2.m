clear all
close all

%% define files
files = {'oct16.mat'; 'feb11.mat'; 'mar15.mat'};
%number of files
nfiles = size(files,1);
%colors for plotting
colors = eye(3);
%number of bootstrap samples when we estimate parameters
nboot = 1000;

%% load data from calibration curve
% I took data from the calibration curve and saved it in curve_data.mat
load curve_data.mat
%x-y parameterizes mean curve
x = interp1(xaxis,time,x);
y = interp1(yaxis,furrow_len,y);
%xmin-ymin parameterizes minimum of curves
xmin = interp1(xaxis,time,xmin);
ymin = interp1(yaxis,furrow_len,ymin);
%xmax-ymax parameterizes maximum of curves
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
%calculate PCs
[I, coeff, V_PCA, D_PCA] = unscramble_pca(dpERK_raw, true);
%calculate mean of data
data_mean = mean(dpERK_raw);
%calculate denoised data
dpERK_denoise = coeff(:,1) * V_PCA(:,1)' + ones(size(coeff,1),1) * data_mean;

%plot data; sort by PCA
figure;
imagesc(dpERK_raw(I,:))
xlabel('position')
ylabel('ordered by PCA')

%plot denoised data; sort by membrane length
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

%% plot PCA modes
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

%% plot eigenvalue spectrum
figure;
plot(diag(D_PCA)/sum(diag(D_PCA)),'.')
xlabel('principal component index')
ylabel('fraction of variance captured')
set(gca,'xlim',[0 10])

%% plot mode with error bars
%compute bootstrap samples and bootstrap confidence intervals
[bootstat_ci, bootstat] = bootci(nboot, {@bootfun_PCA, dpERK_raw});

figure;
errorbar(1:size(bootstat,2),mean(bootstat),mean(bootstat) - bootstat_ci(1,:), bootstat_ci(2,:) - mean(bootstat))
xlabel('position')
ylabel('first PCA mode')

%% compare first PC for each data set; compute error bars using bootstrap
figure;
for i=1:nfiles
    load(files{i,1})    
    [bootstat_ci, bootstat] = bootci(nboot, {@bootfun_PCA, dpERK_raw});
    errorbar(i:3:size(bootstat,2),mean(bootstat(:,i:3:end)),mean(bootstat(:,i:3:end)) - bootstat_ci(1,i:3:end), bootstat_ci(2,i:3:end) - mean(bootstat(:,i:3:end)),'color',colors(i,:))
    hold on
end
xlabel('position')
ylabel('first PCA mode')
legend('data1','data2','data3','location','best')

%% make master mode
%store dpERK profiles from all data sets
total_data = [];
%store times from all data sets
time = [];
%store data set index for each data point
data_ind = [];
for i=1:nfiles
    load(files{i,1})
    
    %add in new data; subtract out the minimum mean profile to "normalize"
    %between data sets
    total_data = [total_data; dpERK_raw-min(mean(dpERK_raw,2))];
    
    %add in times for new data by interpolating using calibration curve
    %data
    time = [time; interp1(y, x, L(:,1), [], 'extrap')];
    
    %store data set index for each point
    data_ind = [data_ind; i*ones(size(dpERK_raw,1),1)];
end
%compute bootstrap samples and bootstrap confidence intervals for all data
[bootstat_ci, bootstat] = bootci(nboot, {@bootfun_PCA, total_data});
master_V = bootfun_PCA(total_data)';

%plot master mode
figure;
errorbar(1:size(bootstat,2),master_V,master_V - bootstat_ci(1,:), bootstat_ci(2,:) -master_V)
xlabel('position')
ylabel('first PCA mode')

%% fit function for PCA coeff vs. time using master mode
%calculate projection coefficients for data
coeff0 = total_data * master_V';
%fit coeff=a*time^b using bootstrap
%power_law_ci stores 95% confidence intervals for power law parameters
%power_law_means stores optimal parameters
[power_law_ci, power_law_data] = bootci(nboot,@power_law_params,time, coeff0);
power_law_means = power_law_params(time, coeff0);

%plot data and fit curve
figure;
scatter(time, coeff0, 50, data_ind, '.')
colormap(colors)
hold on
plot(sort(time), power_law(power_law_means, sort(time)),'k')

xlabel('time from membrane length')
ylabel('first PCA projection coefficient')

%% order data using PCA coefficients and power law function
%plot each data set seperately
for i=1:nfiles  
    %find relevant data
    ind = find(data_ind == i);
    temp_data = total_data(ind,:);
    temp_time = time(ind,:);
    
    %calculate time from power law fit
    pred_time = inv_power_law(power_law_means, coeff0(ind));
    
    %plot data where rows are spaced according to time
    time_range = min(temp_time):0.01:max(temp_time);
    image_data = zeros(length(time_range), size(temp_data,2));
    for j=1:length(time_range)
        [~, idx] = min((pred_time-time_range(j)).^2);
        image_data(j,:) = temp_data(idx,:);
    end
    
    figure;
    imagesc(1:size(image_data,1),time_range,image_data)
    xlabel('position')
    ylabel('time')
end

%% make master mode using only first two data sets; use it to predict third data set
ind = find(data_ind ~= nfiles);
%compute bootstrap samples and bootstrap confidence intervals for all data
[bootstat_ci, bootstat] = bootci(nboot, {@bootfun_PCA, total_data(ind,:)});

%calculate master mode from all data (mean of bootstrap samples)
master_V = mean(bootstat);

%plot master mode
figure;
errorbar(1:size(bootstat,2),master_V,master_V - bootstat_ci(1,:), bootstat_ci(2,:) -master_V)
xlabel('position')
ylabel('first PCA mode')

%fit coeff=a*time^b using bootstrap
[power_law_ci, power_law_data] = bootci(nboot,@power_law_params,time(ind), coeff0(ind));
power_law_means = power_law_params(time(ind), coeff0(ind));

%plot data and fit curve
figure;
scatter(time(ind), coeff0(ind), 50, data_ind(ind), '.')
colormap(colors)
hold on
plot(sort(time), power_law(power_law_means, sort(time)),'k')

xlabel('time from membrane length')
ylabel('first PCA projection coefficient')

%find relevant data
ind2 = find(data_ind == nfiles);
temp_data = total_data(ind2,:);
temp_time = time(ind2,:);

%calculate time from power law fit
pred_time = inv_power_law(power_law_means, coeff0(ind2));

%plot data where rows are spaced according to time
time_range = min(temp_time):0.01:max(temp_time);
image_data = zeros(length(time_range), size(temp_data,2));
for j=1:length(time_range)
    [~, idx] = min((pred_time-time_range(j)).^2);
    image_data(j,:) = temp_data(idx,:);
end

figure;
imagesc(1:size(image_data,1),time_range,image_data)
xlabel('position')
ylabel('time')
