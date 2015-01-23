function [R, W] = compute_pairwise_alignments(images, nrot)
% [R, W] = compute_pairwise_alignments(images, nrot, shift_max, nshifts, display_waitbar)
% compute the pairwise alignments and distances for a set of
% images, to be used in a vector diffusion maps calculation
% images is a set of images, where images(:,:,i) (for grayscale) or
% images(:,:,:,i) (for RGB images) contains the i^th image in the dataset
% nrot is the number of angular discretications to use to compute the
% rotational pairwise alignments; larger values of nrot will increase the
% required computational time, but could increase the accuracy
% shift_max is the largest shift value to consider in the alignments;
% shift_max=0 corresponds to only rotations
% nshifts is the number of shift discretizations to use to compute the
% pairwise alignments; if shift_max=0, the value of nshifts is ignored
%
% R is a matrix containing the pairwise rotation matrices
% W is a matrix containing the minimum pairwise distances

m = size(images, ndims(images));
npixels = size(images, 1);

dim = 2;

% angles to iterate over
theta_vec = linspace(0, 360, nrot+1);
theta_vec = theta_vec(1:end-1);

% store optimal rotations
thetas = zeros(m);
dx = 0;
dy = 0;

R = zeros(dim*m);
W = inf(m);

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);

h_waitbar = waitbar(0,'Computing pairwise alignments...');

for i=1:nrot
    %     waitbar(i/nrot, h_waitbar);
    theta = theta_vec(i);
    
    waitbar(i/nrot, h_waitbar);
    
    % calculate affine transformation
    A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; dx dy 1]);
    
    % apply transform (rotation+translation) images
    images_transformed = imwarp(images, RI, A, 'outputview', RI);
    
    % compute distances between transformed images and initial
    % images
    dist_tmp = pdist2(reshape(double(images), [], m)', reshape(double(images_transformed), [], m)');
    
    % find images where this affine transformation is better than
    % all previous transforms
    idx = find(dist_tmp < W);
    
    % store new distances
    W(idx) = dist_tmp(idx);
    
    % store new rotations and translations
    thetas(idx) = theta;
end

% compute rotation matrices from rotations+translations
for i=1:m
    for j=1:m
        R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix_2d(thetas(i, j));
    end
end

% symmetrize matrices
W = 0.5*(W + W');
R = 0.5*(R + R');

close(h_waitbar);

% compute two-dimensional rotation matrix corresponding to a specific
% rotation
% dtheta is the rotation angle
% R is the rotation matrix
function R = calc_rot_matrix_2d(dtheta)

R = [cosd(dtheta) -sind(dtheta);
    sind(dtheta) cosd(dtheta)];






