clear all
% close all

%%

channel = 1;
npixels = 100;

dt = 0.5;

movies = {'../bomyi_emb01_gast01.avi';...
    '../bomyi_emb02_gast02.avi';...
    '../14_0623/emb01_hisRFP_gastrulation.avi';...
    '../14_0623/emb02_hisRFP_gastrulation.avi';...
    '../14_0624/emb01_hisRFP_gastrulation.avi';...
    '../14_0624/emb02_hisRFP_gastrulation.avi'; ...
     '../0709_emb02_cell_gast.avi'};

nmovies = length(movies);

theta = [-95;
    -55;
    -80;
    -80;
    -95;
    -120;
    -95];

% movie_start = [1; 1; 1; 1; 1; 1; 96];
% movie_end = [20; 18; 8; 7; 0; 0; 0];

movie_start = [20 25 15 10 15 5 110];
  movie_end = [55 60 48 42 42 40 140];


make_subplot = @(i) subplot(2, 4, i);

images = cell(nmovies, 1);
time = cell(nmovies, 1);
nimages = zeros(nmovies, 1);

nmodes = 5;

%%
figure;
for i=1:nmovies
    [images_tmp, time_tmp] = read_video(movies{i}, npixels);
%     images_tmp = images_tmp(:, :, :, movie_start(i):end-movie_end(i));
%     time_tmp = time_tmp(movie_start(i):end-movie_end(i)) - time_tmp(movie_start(i));
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
    
    make_subplot(i);
    imshow(images{i}(:,:, round(nimages(i)/2)));
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
return
%%

train_data = vertcat(PCA_data{:});
train_times = vertcat(time_adjust{:});

[m, n] = size(train_data);

mean_train_data = mean(train_data);

train_data = train_data - repmat(mean_train_data, m, 1);

[U, S, V_train] = svds(train_data, nmodes);

%[b_train, bint_train] = regress(train_times, [train_data*V_train ones(m, 1)]);
mdl = fitlm(train_data*V_train,train_times,'linear');

save('PCA_train_data.mat', 'mean_train_data', 'V_train', 'mdl');

%%

for i=1:nmovies
    [images_tmp, time_tmp] = read_video(movies{i}, npixels);
    images_tmp = images_tmp(:, :, :, movie_start(i):movie_end(i));
    
    time_est = zeros(size(time_adjust{i}));
    time_est_min = zeros(size(time_adjust{i}));
    time_est_max = zeros(size(time_adjust{i}));
    
    for j=1:nimages(i)
        [time_est(j),time_est_min(j),time_est_max(j)] = predict_time_image(imrotate(images_tmp(:,:,channel,j), theta(i), 'crop'));
    end
    
    figure;
    subplot(1,3,1)
    plot(time_adjust{i}, time_est,'.')
    hold on
    plot([0 15], [0 15],'-r')

    subplot(1,3,2)
    plot(abs(time_adjust{i} - time_est),'.')

    subplot(1,3,3)
     plot(time_est_max - time_est_min)

end

