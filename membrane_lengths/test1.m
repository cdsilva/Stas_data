clear all
close all

imshow(imread('Ingression_kinetics_Anna Sokac.tif'))

load curve_data.mat
x = interp1(xaxis,time,x);
y = interp1(yaxis,furrow_len,y);
xmin = interp1(xaxis,time,xmin);
ymin = interp1(yaxis,furrow_len,ymin);
xmax = interp1(xaxis,time,xmax);
ymax = interp1(yaxis,furrow_len,ymax);

%lerr = y-interp1(xmin, ymin, x, [], 'extrap');
%uerr = interp1(xmax, ymax, x, [], 'extrap')-y;

figure;
%errorbar(x,y,y-interp1(xmin, ymin, x),interp1(xmax,ymax,x)-y)
%errorbar(x,y,lerr,uerr)
%hold on
%plot(xmin,ymin)
%plot(xmax,ymax)

%figure;
plot(x,y)
hold on
plot(xmin, ymin)
plot(xmax, ymax)
set(gca,'ydir','reverse')

files = {'oct16.mat'; 'feb11.mat'; 'mar15.mat'};
nfiles = size(files,1);

mfig = nfiles;
nfig = 5;

figure;
for i=1:nfiles
    load(files{i,1})
    
    subplot(mfig, nfig, 1+(i-1)*nfig)
    imagesc(dpERK_sort)
    xlabel('position')
    ylabel('membrane order')
    title('sorted by furrow length')
    
    [I, coeff, V_PCA, D_PCA] = unscramble_pca(dpERK_raw);
    
    subplot(mfig, nfig, 2+(i-1)*nfig)
    imagesc(dpERK_raw(I,:))
    xlabel('position')
    ylabel('PCA order')
    title('sorted by PCA')
    
    subplot(mfig, nfig, 3+(i-1)*nfig)
    imagesc((dpERK_raw(I,:)-dpERK_sort).^2)
    xlabel('position')
    ylabel('time')
    title('squared error between orderings')
    
    time = interp1(y, x, L(:,1), [], 'extrap');
    time_max = interp1(ymin, xmin, L(:,1), [], 'extrap');
    time_min = interp1(ymax, xmax, L(:,1), [], 'extrap');
    
    traj = zeros(length(time),1);
    traj(1) = time(I(1)) - time_min(I(1));
    for j=2:length(time)
        traj(j) = max(traj(j-1),time_min(I(j)));
    end
    
    subplot(mfig, nfig, 4+(i-1)*nfig);
    plot(time(I),'.')
    hold on
    plot(time_min(I))
    plot(time_max(I))
    hold on
    %plot(traj,'r')
    %scatter(1:length(time), time(I), 100, traj > time_max(I)+time(I), '.')
    %colormap(winter)
    xlabel('PCA order')
    ylabel('time from furrow length')
    axis([-inf inf -inf inf])
    %view(90,-90)
    %set(gca,'xdir','reverse')
    
%     subplot(mfig, nfig, 5+(i-1)*nfig);
%     plot(time, coeff(:,1),'.')

    disp(files{i})
    disp(num2str(D_PCA(1,1)));
    disp(num2str(D_PCA(2,2)));
    disp(num2str(D_PCA(1,1)/sum(diag(D_PCA))));
    
    
    subplot(mfig, nfig, 5+(i-1)*nfig);
    [~, rank_membrane] = sort(L(:,2));
    [~, rank_PCA] = sort(I);
    plot(rank_membrane, rank_PCA,'.')
    disp(num2str(corr(L(:,1),coeff(:,1),'Type','Kendall')))
    

    

end