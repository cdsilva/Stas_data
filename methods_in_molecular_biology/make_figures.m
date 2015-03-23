% clear all
close all

%%

channel = 1;
npixels = 100;

dt = 0.5;

dir = '../live_imaging';
movies = {'bomyi_emb01_gast01.avi';...
    'bomyi_emb02_gast02.avi';...
    '14_0623/emb01_hisRFP_gastrulation.avi';...
    '14_0623/emb02_hisRFP_gastrulation.avi';...
    '14_0624/emb01_hisRFP_gastrulation.avi';...
    '14_0624/emb02_hisRFP_gastrulation.avi'; ...
    '0709_emb02_cell_gast.avi'};

nmovies = length(movies);

theta = [-95;
    -55;
    -80;
    -80;
    -95;
    -120;
    -95];

movie_start = [20 25 15 10 15 5 110];
movie_end = [55 60 48 42 42 40 140];

images = cell(nmovies, 1);
time = cell(nmovies, 1);
nimages = zeros(nmovies, 1);

nmodes = 5;

%%
for i=1:nmovies
    [images_tmp, time_tmp] = read_video(sprintf('%s/%s', dir, movies{i}), npixels);
    images_tmp = images_tmp(:, :, :, movie_start(i):movie_end(i));
    time_tmp = time_tmp(movie_start(i):movie_end(i)) - time_tmp(movie_start(i));

    time{i} = time_tmp * dt;

    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);

    nimages(i) = length(time{i});

    for j=1:nimages(i)
        images_tmp(:, :, j) = image_normalize(images_tmp(:, :, j), theta(i));
    end

    images{i} = images_tmp;

end

%%

PCA_data = cell(nmovies, 1);
for i=1:nmovies
    PCA_data{i} = create_PCA_data(images{i});
end

%%

shift_times = estimate_shift_times(PCA_data, time, nmodes);

time_adjust = time;
for i=1:nmovies
    time_adjust{i} = time_adjust{i} - shift_times(i);
end

%%

train_data = vertcat(PCA_data{:});
train_times = vertcat(time_adjust{:});

[m, n] = size(train_data);

mean_train_data = mean(train_data);

train_data = train_data - repmat(mean_train_data, m, 1);

[U, S, V_train] = svds(train_data, nmodes);


%% store eigenimages

eigenimages = zeros(npixels, npixels, nmodes);

for i=1:nmodes
    if mean(V_train(:,i)) < 0
        V_train(:,i) = -V_train(:,i);
    end
    eigenimages(:,:,i) = reshape(V_train(:,i), [npixels npixels]);
    %     figure;
    %     imagesc(eigenimages(:,:,i))
    %     colormap(gray)
    %     axis off
    %     axis square
    %
end

return
%% regression with one variable

i = 1;

y = train_times - min(train_times);
x = train_data*V_train(:,1:i);

mdl = fitlm(x, y,'linear');
[Y,~] = predict(mdl ,sort(x));
[pred_time,~] = predict(mdl ,x);

make_fig(6,6);
scatter(x(:,1),y,50, 'b', '.')
hold on
axis square
plot(sort(x), Y, '-r')
xlabel('variable 1')
ylabel('time')
print(gcf, '-dpdf', 'figures/corr1')

make_fig(8,8);
plot(y, pred_time, '.')
xlabel('time')
ylabel('predicted time')
axis equal
axis square
print(gcf, '-dpdf', 'figures/errors1')

%% regression with two variables

i = 2;

x = train_data*V_train(:,1:i);

make_fig(10, 8);
scatter(x(:,1),x(:,2),50, y, '.')
xlabel('variable 1')
ylabel('variable 2')
axis tight
axis square
h = colorbar;
set(get(h,'xlabel'),'String', 'time');
print(gcf, '-dpdf', 'figures/corr2_data')

x1 = linspace(min(x(:,1)), max(x(:,1)), 100);
x2 = linspace(min(x(:,2)), max(x(:,2)), 100);
[X1, X2] = meshgrid(x1, x2);

mdl = fitlm(x, y,'linear');
[Y,~] = predict(mdl ,[X1(:) X2(:)]);
[pred_time,~] = predict(mdl ,x);

make_fig(6,6);
% scatter(X1(:), X2(:), 50, Y, '.')
imagesc(x1, x2, reshape(Y, size(X1)))
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
% hold on
% scatter(x(:,1),x(:,2),50, y, '.')
% plot(x(:,1),x(:,2),'ow','markersize',2,'linewidth', 0.2)
xlabel('variable 1')
ylabel('variable 2')
axis square
print(gcf, '-dpdf', 'figures/corr2_fit')

make_fig(6,6);
plot(y, pred_time, '.')
axis equal
axis square
xlabel('time')
ylabel('predicted time')
print(gcf, '-dpdf', 'figures/errors2')

return

%% regression with three-five variables

for i=3:5
    
    x = train_data*V_train(:,1:i);
    
    mdl = fitlm(x, y,'linear');
    [pred_time,~] = predict(mdl ,x);
    
    make_fig(5,5);
    plot(y, pred_time, '.')
    xlabel('time')
    ylabel('predicted time')
    print(gcf, '-dpdf', sprintf('figures/errors%d', i))
    
end


for j=1:5
    image_scales = [min(x(:,j)) mean(x(:,j)) max(x(:,j))];
    
    images_to_plot = repmat(eigenimages(:,:,j), [1 1 length(image_scales)]);
    for i=1:length(image_scales)
        images_to_plot(:,:,i) = images_to_plot(:,:,i) * image_scales(i);
    end
    clim = [min(images_to_plot(:)) max(images_to_plot(:))];
    
    for i=1:length(image_scales)
        make_fig(2,2);
        imagesc(images_to_plot(:,:,i), clim)
        colormap(gray)
        
        axis off
%         axis square
        set(gca,'position',[0 0 1 1],'units','normalized')

        print(gcf, '-dpdf', sprintf('figures/eigenimage%d_%d', j, i))
    end
end

