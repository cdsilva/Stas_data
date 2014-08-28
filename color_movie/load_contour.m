function [images, time] = load_contour(npixels)
% movie with countours

nspline_points = 3000;

fid = fopen('spline_x.txt');
A = textscan(fid, '%f');
fclose(fid);
x = reshape(A{1}, nspline_points+1, []);
x = x(2:end, :);

fid = fopen('spline_y.txt');
A = textscan(fid, '%f');
fclose(fid);
y = reshape(A{1}, nspline_points+1, []);
y = y(2:end, :);

ncontour_images = size(y, 2);

images = zeros(npixels, npixels, ncontour_images, 'uint8');

h = figure;
for i=1:ncontour_images
    plot(x(:,i), y(:,i), 'linewidth', 20, 'color','k')
    axis([0 1000 0 1000])
    axis square
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    A = getframe;
    images(:,:,i) = imresize(imcomplement(rgb2gray(A.cdata(2:end-1,2:end-1,:))), [npixels npixels]);
    clf
end
close(h);

movie_angle = -125;
images = imrotate(images, movie_angle, 'crop');

dt = 0.5;
time = dt * (1:ncontour_images)';

for i=1:ncontour_images
    
    im_tmp = images(:,:,i);
    
    row_sum = sum(im_tmp);
    idx1 = find(row_sum > 0, 1, 'first');
    idx2 = find(row_sum > 0, 1, 'last');
    
    col_sum = sum(im_tmp, 2);
    idx3 = find(col_sum > 0, 1, 'first');
    idx4 = find(col_sum > 0, 1, 'last');
    
    images(:,:,i) = imresize(im_tmp(idx3:idx4, idx1:idx2), [npixels npixels]);
    
end



