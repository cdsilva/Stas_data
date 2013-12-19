function [R, W, angles] = align_data_nosph(images, angle_proj)

dim = 3;

m = size(images, 3);
n = size(images, 1);

R = zeros(dim*m);
angles = zeros(m, m, 3);
W = inf(m);

nrot = 100;
theta = linspace(0, 2*pi, nrot+1);
theta = theta(1:end-1);

nshift = 17;
shift_max = 8;
if mod(nshift, 2) ~= 1 || mod(shift_max, (nshift-1)/2) ~= 0 
    disp('nshift is not valid')
    return
end
shifts = linspace(-shift_max, shift_max, nshift);

for i=1:m
    imagei = uint8(images(:,:,i));
    for k=1:nrot
        imagei_tmp = imrotate(imagei, (180/pi)*theta(k), 'crop');
        for k1=1:nshift
            for k2=1:nshift
                imagei_tmp2 = circshift(imagei_tmp, [shifts(k1) shifts(k2)]);
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
        %if W(i,j) < W(j,i)
            W(j,i) = W(i,j);
            alpha = angles(i, j, 1);
            beta = angles(i, j, 2) / n * angle_proj;
            gamma = angles(i, j, 3) / n * angle_proj;
%         else
%             W(i,j) = W(j,i);
%             alpha = angles(j, i, 1);
%             beta = angles(j, i, 2) / n * pi/4;
%             gamma = angles(j, i, 3) / n * pi/4;
%         end
        %R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = rot_matrix([alpha, beta, gamma]);
        R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = rot_matrix2(alpha, beta, gamma);
        R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i) = R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j)';
        
    end
end

function R = rot_matrix2(alpha, beta, gamma)

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



