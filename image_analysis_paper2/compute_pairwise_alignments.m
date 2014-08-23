function [R, W] = compute_pairwise_alignments(images, nrot, nshifts, shift_max)

m = size(images, ndims(images));
npixels = size(images, 1);

% only rotations
if nshifts > 1
    dim = 3;
    angle_proj = 30;
    shifts = linspace(-shift_max, shift_max, nshifts);
    % rotations + translations
else
    dim = 2;
    nshifts = 1;
    shifts = 0;
end

theta_vec = linspace(0, 360, nrot+1);
theta_vec = theta_vec(1:end-1);

thetas = zeros(m);
dx = zeros(m);
dy = zeros(m);

R = zeros(dim*m);
W = inf(m);

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);

h_waitbar = waitbar(0,'Computing pairwise alignments...');
for i=1:nrot
%     waitbar(i/nrot, h_waitbar);
    theta = theta_vec(i);
    for j=1:nshifts
        for k=1:nshifts
            waitbar(((i-1)*nshifts^2+(j-1)*nshifts+k)/(nrot*nshifts^2), h_waitbar);
            A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; shifts(j) shifts(k) 1]);
            images_transformed = imwarp(images, RI, A, 'outputview', RI);
            dist_tmp = pdist2(reshape(double(images), [], m)', reshape(double(images_transformed), [], m)');
            idx = find(dist_tmp < W);
            W(idx) = dist_tmp(idx);
            thetas(idx) = theta;
            dx(idx) = shifts(j);
            dy(idx) = shifts(k);
        end
    end
end

if dim == 3
    for i=1:m
        for j=1:m
            R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix_3d(thetas(i, j), dx(i, j), dy(i, j), angle_proj);
        end
    end
else %dim == 2
    for i=1:m
        for j=1:m
            R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix_2d(thetas(i, j));
        end
    end
end

W = 0.5*(W + W');
R = 0.5*(R + R');

close(h_waitbar);

function R = calc_rot_matrix_3d(dtheta, dx, dy, angle_proj)

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

function R = calc_rot_matrix_2d(dtheta)

R = [cosd(dtheta) -sind(dtheta);
    sind(dtheta) cosd(dtheta)];






