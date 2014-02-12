function [R, W, angles] = align_data_nosph2(images, angle_proj, shift_max, delta_shift)

dim = 3;

m = size(images, 3);
n = size(images, 1);

R = zeros(dim*m);
angles = zeros(m, m, 3);
W = inf(m);

nrot = 20;
theta = linspace(0, 2*pi, nrot+1);
theta = theta(1:end-1);

shifts = (-shift_max:delta_shift:shift_max) * angle_proj / n;
nshift = length(shifts);

all_angles = zeros(nrot*nshift*nshift, 3);
all_rot_images = zeros(nrot*nshift*nshift, n*n, m);

for i=1:m
    imagei = uint8(images(:,:,i));
    angle_idx = 1;
    for k=1:nrot
        imagei_tmp = rotate_image(imagei, theta(k));
        for k1=1:nshift
            for k2=1:nshift
                imagei_tmp2 = translate_image(imagei_tmp, shifts(k1), shifts(k2), angle_proj);
                all_rot_images(angle_idx, :, i) = reshape(double(imagei_tmp2), 1, []);
                all_angles(angle_idx, 1) = theta(k);
                all_angles(angle_idx, 2) = shifts(k1);
                all_angles(angle_idx, 3) = shifts(k2);
                angle_idx = angle_idx + 1;
            end
        end
    end
end


for i=1:m
    for j=1:i-1
        [d2, idx] = min(sum((all_rot_images(:,:,i) - repmat(reshape(double(images(:,:,j)), 1, []), nrot*nshift*nshift, 1)).^2, 2));
        W(i,j) = d2;
        angles(i, j, :) = all_angles(idx, :);
    end
end


for i=1:m
    for j=1:i-1
        W(j,i) = W(i,j);
        R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix(angles(i, j, 1), angles(i, j, 2), angles(i, j, 3));
        R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i) = R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j)';
    end
end





