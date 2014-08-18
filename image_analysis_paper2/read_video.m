function [images, times] = read_video(avi_file, npixels)

vidobj = VideoReader(avi_file);

num_images = vidobj.NumberOfFrames;

A = read(vidobj);
images = imresize(A, [npixels npixels]);

times = (1:num_images)';

