function [R, W] = compute_pairwise_alignments_color(images, angle_proj, shift_max, shift_step)

dim = 3;

m = size(images, 4);
npixels = size(images, 1);

R = zeros(dim*m);
W = inf(m);

nrot = 20;
theta_vec = linspace(0, 360, nrot+1);
theta_vec = theta_vec(1:end-1);

shifts = -shift_max:shift_step:shift_max;
nshift = length(shifts);

thetas = zeros(m);
dx = zeros(m);
dy = zeros(m);

parfor i=1:m
    i
    imagei = images(:,:,:,i);
    W_tmp = inf(1, m);
    thetas_tmp = zeros(m, 1);
    dx_tmp = zeros(m, 1);
    dy_tmp = zeros(m, 1);
    for k=1:nrot
        imagei_tmp = imrotate(imagei, theta_vec(k), 'crop');
        for k1=1:nshift
            for k2=1:nshift
                imagei_tmp2 = circshift(imagei_tmp, [shifts(k1) shifts(k2) 0]);
                if shifts(k1) > 0
                    imagei_tmp2(1:shifts(k1), :, :) = 0;
                end
                if shifts(k1) < 0
                    imagei_tmp2(end+shifts(k1)+1:end, :, :) = 0;
                end
                if shifts(k2) > 0
                    imagei_tmp2(:, 1:shifts(k2), :) = 0;
                end
                if shifts(k2) < 0
                    imagei_tmp2(:,end+shifts(k2)+1:end, :) = 0;
                end
                
                for j=1:i-1
                    imagej = images(:,:,:,j);
                    d2 = sum((double(imagei_tmp2(:))-double(imagej(:))).^2);
                    %                     if d2 < W(i, j)
                    %                         W(i,j) = d2;
                    %                         thetas(j, i) = theta_vec(k);
                    %                         dx(j, i) = shifts(k1);
                    %                         dy(j, i) = shifts(k2);
                    %                     end
                    if d2 < W_tmp(j)
                        W_tmp(j) = d2;
                        thetas_tmp(j) = theta_vec(k);
                        dx_tmp(j) = shifts(k1);
                        dy_tmp(j) = shifts(k2);
                    end
                end
            end
        end
    end
    W(i, :) = W_tmp;
    thetas(:, i) = thetas_tmp;
    dx(:, i) = dx_tmp;
    dy(:, i) = dy_tmp;
end

for i=1:m
    for j=1:i-1
        W(j,i) = W(i,j);
        R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i) = calc_rot_matrix(thetas(j, i), dx(j, i), dy(j, i), npixels, angle_proj);
        R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i)';
    end
end

function R = calc_rot_matrix(dtheta, dx, dy, npixels, angle_proj)

alpha = dtheta * pi / 180;
beta = dx /npixels * angle_proj;
gamma = dy / npixels * angle_proj;

Rx = [1 0 0;
    0 cos(alpha) -sin(alpha);
    0 sin(alpha) cos(alpha)];

Ry = [cos(beta) 0 sin(beta);
    0 1 0;
    -sin(beta) 0 cos(beta)];

Rz = [cos(gamma) -sin(gamma) 0;
    sin(gamma) cos(gamma) 0;
    0  0 1];

R = Rz * Ry * Rx;





