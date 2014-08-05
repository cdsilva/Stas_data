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

for i=1:nmovies
    figure;
    subplot_dim1 = ceil(sqrt(nimages(i)));
    subplot_dim2 = ceil(nimages(i) / subplot_dim1);
    
    for j=1:nimages(i)
        subplot(subplot_dim1, subplot_dim2, j)
        imshow(images{i}(:,:,j));
    end
    
end

%%

PCA_data = cell(nmovies, 1);
for i=1:nmovies
    PCA_data{i} = create_PCA_data(images{i});
end



%%

nmodes = 8;

shift_times = estimate_shift_times(PCA_data, time, nmodes);

time_adjust = time;
for i=1:nmovies
    time_adjust{i} = time_adjust{i} - shift_times(i);
end

%%
bias = zeros(nmovies);
variance = zeros(nmovies);

for i=1:nmovies
    [B, FitInfo] = lasso(PCA_data{i},time_adjust{i}, 'lambda', 0.5);
    
    figure;
    for j=1:nmovies
        make_subplot(j)
        plot(time_adjust{j}, PCA_data{j}*B+FitInfo.Intercept, '.')
        hold on
        plot(time_adjust{j}, time_adjust{j}, '-r');
        
        bias(i, j) = mean(time_adjust{j} - (PCA_data{j}*B+FitInfo.Intercept));
        variance(i, j) = mean((time_adjust{j} - (PCA_data{j}*B+FitInfo.Intercept)).^2);
    end
end

figure;
imagesc(bias)
colorbar

figure;
imagesc(sqrt(variance))
colorbar

figure; 
imagesc(reshape(B~=0, [100 100]))
