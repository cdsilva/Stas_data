function [images, time] = load_movie(npixels)
% movie with nuclei (red)

movie_file = '../live_imaging/14_0624/emb02_cell_gast.avi';
movie_angle = -120;
movie_end = 142;

channel = 1;

dt = 0.5;

[images, time] = read_video(movie_file, npixels);
images = images(:,:,:,1:movie_end);
time = dt * time(1:movie_end);
nmovie_images = length(time);

images = imrotate(images, movie_angle, 'crop');
images(:,:,2,:) = 0;

for i=1:nmovie_images
    images(:,:,:,i) = crop_hires_image(images(:,:,:,i), channel);
end