clear all
close all

%% Basic Algorithm
% The basic steps are as follows:
%
% # Read in images 
% # Preprocess the images (smooth and/or equalize channels, mean-center, etc.)
% # Run vector diffusion maps to register and order preprocessed images

%% Read Images
% The first step is to read in images.
% This is done using the |read_images| function 
% For two-dimensional images, it requires that all the images are stored in
% a single directory, all with a consistent prefix, and indexed with two
% digits starting from |01|. 
% For three-dimensional z-stack, it assumes each z-stack is stored in a
% single directory. 
% All of the z-stack directories are required to begin with a consistent prefix, and be indexed with two
% digits starting from |01|. 
% Similarly, all of the images are required to begin with a consistent prefix, and be indexed with two
% digits starting from |01| for each z-stack.

% directory where images are stored
image_dir = 'drosophila_fixed_images';

% prefix for each image
image_name = 'emb';

% image type/extension
image_ext = 'tif';

% no stack name or number of iamges in a stack are required for 2D images
stack_name = '';
nstack = 0;

% number of images
nimages = 120;

% dimension of images (dim=2 indicates standard 2D images, rather than
% z-stacks)
dim = 2;

% read in images
% images are stored in the variable images_raw
[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, ...
    stack_name, nimages, nstack, dim);

%% Show images
% An optional step is to now show the images using the |plot_images|
% function.

% plot the images
% images_raw are the images to plot 
% (returned from the read_images function)
% dim is the dimension of the images 
% (dim=2 indicates standard 2D images, rather than z-stacks)
plot_images(images_raw, dim)

%% Preprocess Images
% Now, we must preprocess the images before registration and ordering, to
% remove any imaging and/or experimental artifacts.
% We do this via the |apply_image_functions| function
% For this particular imaging data set, there are three channels. 
% DAPI, a nuclear stain, is in the first (red) channel.
% dpERK is in the second (green) channel.
% Dl is in the third (blue) channel.

% number of pixels (fewer pixels means a shorter runtime)
npixels = 100;

% channel weights
% because the nuclear signal is wider spread throughout the image, as well
% as noiser, we choose to scale the first channel by half
channel_weight = [0.5 1 1];

% we blur the image slightly (5%) so that small-scale effects of individual
% nuclei are smoothened
channel_blur = [0.05 0.05 0.05];

% because the absolute intensity of the nuclear signal is not informative
% or important (only the presence or absence of signal is important, 
% indicating the presence or absence of nuclei), we normalize the first channel
% We do not normalize the other two channels, as changes in intensity are
% informative to the developmental processes.
channel_normalize = [1 0 0];

% we use the nuclear channel to mean center the images, as this channel
% parameterizes/captures the entire embryo
channel_mean_center = [1 0 0];

% we resize the images to all be (approximately) the same size, to remove
% any variations due to size effects, as size effects are known to not be
% important for this developmental system
resize_image = true;

% we then apply these image functions of normalization, blurring,
% reweighting, and mean-centering
images = apply_image_functions(images_raw, npixels, dim, channel_weight, ...
    channel_blur, channel_normalize, channel_mean_center, resize_image);

% plot the images (optional)
plot_images(images, dim)

%% Calculate pairwise alignments
% We now need to calculate the angles needed to align _pairs_ of images

% angular discretization when computing pairwise aligments
% this means we search for pairwise aligmemnts over 10 degree increments
ang_dis = 10;

% compute the pairwise alignments 
% images are the preprocessed images
% ang_dis is the angular discretization
% R and W store the pairwise alignments and distances, respectively, for
% vector diffusion maps
[R, W] = compute_pairwise_alignments(images, ang_dis);

%% Apply vector diffusion maps
% We can now use vector diffusion maps to register and order the images. 

% ncomps is the number of components to compute
% we only compute 1 coordinate because we only need to order the images
% (i.e., sort by the first coordinate)
ncomps = 1;

% epsilon scale for diffusion maps kernel
% eps_scale = 0.25 means that the epsilon in the diffusion maps kernel is
% 1/4 of the median of the pairwise distances between data points
eps_scale = 0.25;

% vector diffusion maps calculates optimal rotations and embedding
% coordinate
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images
images_registered = register_all_images(images, R_opt);

% order registered images by embedding coordinate
images_analyzed = order_all_images(images_registered, embed_coord);

% plot the images (optional)
plot_images(images_analyzed, dim)

