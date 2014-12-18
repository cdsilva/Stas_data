function [R, W] = compute_pairwise_alignments(images, nrot, shift_max, nshifts, display_waitbar)
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
% display_waitbar indicates whether to display a waitbar indicating the
% progress of the computation; display_waitbar=true is the default, if no
% value is supplied
%
% R is a matrix containing the pairwise rotation matrices
% W is a matrix containing the minimum pairwise distances

if nargin < 5
    display_waitbar = true;
end

m = size(images, ndims(images));
npixels = size(images, 1);

% only rotations
if shift_max > 0
    dim = 3;
    shifts = linspace(-shift_max, shift_max, nshifts);
% rotations + translations
else
    dim = 2;
    nshifts = 1;
    shifts = 0;
end

% angles to iterate over
theta_vec = linspace(0, 360, nrot+1);
theta_vec = theta_vec(1:end-1);

% store optimal translations and rotations
thetas = zeros(m);
dx = zeros(m);
dy = zeros(m);

R = zeros(dim*m);
W = inf(m);

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);

if display_waitbar
    h_waitbar = waitbar(0,'Computing pairwise alignments...');
end
for i=1:nrot
%     waitbar(i/nrot, h_waitbar);
    theta = theta_vec(i);
    for j=1:nshifts
        for k=1:nshifts
            
            if display_waitbar
                waitbar(((i-1)*nshifts^2+(j-1)*nshifts+k)/(nrot*nshifts^2), h_waitbar);
            end
            
            % calculate affine transformation
            A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; shifts(j) shifts(k) 1]);
            
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
            dx(idx) = shifts(j);
            dy(idx) = shifts(k);
        end
    end
end

% compute rotation matrices from rotations+translations
if dim == 3
    for i=1:m
        for j=1:m
            R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix_3d(thetas(i, j), dx(i, j), dy(i, j));
        end
    end
else %dim == 2
    for i=1:m
        for j=1:m
            R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix_2d(thetas(i, j));
        end
    end
end

% symmetrize matrices
W = 0.5*(W + W');
R = 0.5*(R + R');

if display_waitbar
    close(h_waitbar);
end

% compute three-dimensional rotation matrix corresponding to a specific
% translation and rotation
% dtheta is the rotation angle
% dx is the x-translation (in fraction of the image)
% dy is the y-translation (in fraction of the image)
% R is the rotation matrix
function R = calc_rot_matrix_3d(dtheta, dx, dy)

angle_proj = 20;

alpha = dtheta;
beta = dx * angle_proj;
gamma = dy * angle_proj;

Rx = [1 0 0;
    0 cosd(alpha) -sind(alpha);
    0 sind(alpha) cosd(alpha)];

Ry = [cosd(beta) 0 sind(beta);
    0 1 0;
    -sind(beta) 0 cosd(beta)];

Rz = [cosd(gamma) -sind(gamma) 0;
    sind(gamma) cosd(gamma) 0;
    0  0 1];

R = Rz * Ry * Rx;

% compute two-dimensional rotation matrix corresponding to a specific
% rotation
% dtheta is the rotation angle
% R is the rotation matrix
function R = calc_rot_matrix_2d(dtheta)

R = [cosd(dtheta) -sind(dtheta);
    sind(dtheta) cosd(dtheta)];






