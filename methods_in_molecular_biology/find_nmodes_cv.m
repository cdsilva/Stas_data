clear all
close all

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

movie_start = [20 25 15 10 15 5 110];
  movie_end = [55 60 48 42 42 40 140];

make_subplot = @(i) subplot(2, 4, i);

images = cell(nmovies, 1);
images_raw = cell(nmovies, 1);
time = cell(nmovies, 1);
nimages = zeros(nmovies, 1);

%%
figure;
for i=1:nmovies
    [images_tmp, time_tmp] = read_video(movies{i}, npixels);
    images_tmp = images_tmp(:, :, :, movie_start(i):movie_end(i));
    time_tmp = time_tmp(movie_start(i):movie_end(i)) - time_tmp(movie_start(i));
    
    time{i} = time_tmp * dt;
    
    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);
    images_tmp_raw = zeros(size(images_tmp), 'uint8');
    
    nimages(i) = length(time{i});
    
    for j=1:nimages(i)
        images_tmp(:, :, j) = image_normalize(images_tmp(:, :, j), theta(i));
        images_tmp_raw(:, :, j) = imrotate(images_tmp(:, :, j), theta(i), 'crop');
    end
    
    images{i} = images_tmp;
        images_raw{i} = images_tmp_raw;

    make_subplot(i);
    imshow(images{i}(:,:, round(nimages(i)/2)));
end

%%

PCA_data = cell(nmovies, 1);
for i=1:nmovies
    PCA_data{i} = create_PCA_data(images{i});
end

%%

nmodes_all = 1:10;
cv_err = zeros(size(nmodes_all));

for j = 1:length(nmodes_all);
    j
    
    nmodes = nmodes_all(j);
    
    shift_times = estimate_shift_times(PCA_data, time, nmodes);
    
    time_adjust = time;
    for i=1:nmovies
        time_adjust{i} = time_adjust{i} - shift_times(i);
    end
    
    bias = zeros(nmovies, 1);
    variance = zeros(nmovies, 1);
    for i=1:nmovies
        
        data_tmp = vertcat(PCA_data{1:nmovies ~= i});
        time_tmp = vertcat(time_adjust{1:nmovies ~= i});
        
        pred_time = predict_times_PCA(data_tmp, time_tmp, PCA_data{i}, nmodes);
        
        bias(i) = mean(pred_time - time_adjust{i});
        variance(i) = mean((pred_time - time_adjust{i}).^2);
    end
    
    cv_err(j) = sqrt(mean(variance));
end

figure;
plot(nmodes_all, cv_err, '.')
xlabel('number of PCA modes')
ylabel('leave one out cross-validation error (standard deviation in minutes)')
