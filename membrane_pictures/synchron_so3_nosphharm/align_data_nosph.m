function [R, W, angles] = align_data_nosph(images, angle_proj)

dim = 3;

m = size(images, 3);
n = size(images, 1);

R = zeros(dim*m);
angles = zeros(m, m, 3);
W = inf(m);

nrot = 20;
theta = linspace(0, 2*pi, nrot+1);
theta = theta(1:end-1);
    
nshift = 11;
shift_max = 10;
if mod(nshift, 2) ~= 1 || mod(shift_max, (nshift-1)/2) ~= 0
    disp('nshift is not valid')
    return
end
shifts = linspace(-shift_max, shift_max, nshift) * angle_proj / n;

for i=1:m
    imagei = uint8(images(:,:,i));
    for k=1:nrot
        imagei_tmp = rotate_image(imagei, theta(k));
        for k1=1:nshift
            for k2=1:nshift
                imagei_tmp2 = translate_image(imagei_tmp, shifts(k1), shifts(k2), angle_proj);
                for j=1:i-1
                    d2 = sum(sum((images(:,:,j) - double(imagei_tmp2)).^2));
                    if d2 < W(i, j)
                        W(i,j) = d2;
                        angles(i, j, 1) = theta(k);
                        angles(i, j, 2) = shifts(k1);
                        angles(i, j, 3) = shifts(k2);
                    end
                end
            end
        end
    end
end

for i=1:m
    for j=1:i-1
        W(j,i) = W(i,j);        
        R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = calc_rot_matrix(angles(i, j, 1), angles(i, j, 2), angles(i, j, 3));
        R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i) = R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j)';
    end
end





