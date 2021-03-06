clear all
close all

%%

channel = 1;
npixels = 100;

dt = 0.5;

crop_threshold = 1e3;

% image_equalize = @(image, rot_angle) crop_image(adapthisteq(imrotate(image, rot_angle, 'crop')), crop_threshold, 0);
image_equalize = @(image, rot_angle) crop_image(adapthisteq(imrotate(image, rot_angle, 'crop'), 'cliplimit',0.005), crop_threshold, 0);
image_symmetrize = @(image) imlincomb(0.5, image, 0.5, fliplr(image));
image_fn = @(image, rot_angle) adapthisteq(image_symmetrize(image_equalize(image, rot_angle)));

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

movie_start = [1; 1; 1; 1; 1; 1; 96];
movie_end = [20; 18; 8; 7; 0; 0; 0];

make_subplot = @(i) subplot(2, 4, i);

images = cell(nmovies, 1);
time = cell(nmovies, 1);
nimages = zeros(nmovies, 1);

%%
figure;
for i=1:nmovies
    [images_tmp, time_tmp] = read_video(movies{i}, npixels);
    images_tmp = images_tmp(:, :, :, movie_start(i):end-movie_end(i));
    time_tmp = time_tmp(movie_start(i):end-movie_end(i)) - time_tmp(movie_start(i));
    
    time{i} = time_tmp * dt;
    
    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);
    
    nimages(i) = length(time{i});
    
    for j=1:nimages(i)
        images_tmp(:, :, j) = image_fn(images_tmp(:, :, j), theta(i));
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
