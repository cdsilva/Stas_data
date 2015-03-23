function [images, times] = read_video(avi_file, npixels)

vidobj = VideoReader(avi_file);

num_images = vidobj.NumberOfFrames;

nchannels = 3;
images = zeros(npixels, npixels, nchannels, num_images, 'uint8');

for i = 1:num_images
    A = read(vidobj, i);
    A = imresize(A, [npixels npixels]);
    images(:, :, :, i) = A;
end

times = (1:num_images)';

