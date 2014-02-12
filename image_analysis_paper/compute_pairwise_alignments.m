function [R, W] = compute_pairwise_alignments(images, xmax, ymax, shift_max, shift_step)

nrot = 20;
theta = linspace(0, 360, nrot + 1);
theta = theta(1:end-1);

shifts = -shift_max:shift_step:shift_max;
nshift = length(shifts);

dim = 3;

m = size(images, 3);
npixels = size(images, 1);

R = zeros(dim*m);
W = inf(m);
theta_opt = zeros(m);
xshift_opt = zeros(m);
yshift_opt = zeros(m);


for i=1:m
    imagei = uint8(images(:,:,i));
    for k=1:nrot
        imagei_tmp = imrotate(imagei, theta(k), 'crop');
        for k1=1:nshift
            for k2=1:nshift
                imagei_tmp2 = circshift(imagei_tmp, [shifts(k1) shifts(k2)]);
                for j=1:i-1
                    d2 = sum(sum((images(:,:,j) - double(imagei_tmp2)).^2));
                    if d2 < W(i, j)
                        W(i,j) = d2;
                        theta_opt(i, j) = theta(k);
                        xshift_opt(i, j) = shifts(k1);
                        yshift_opt(i, j) = shifts(k2);
                    end
                end
            end
        end
    end
end

for i=1:m
    for j=1:i-1
        Rtmp = calc_rot_matrix(xmax, ymax, npixels, xshift_opt(i,j), yshift_opt(i,j), theta_opt(i,j));
        R(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j) = Rtmp;
        R(dim*(j-1)+1:dim*j,dim*(i-1)+1:dim*i) = Rtmp';
        W(j,i) = W(i,j);
    end
end

function R = calc_rot_matrix(xmax, ymax, npixels, dxpixels, dypixels, dtheta)

% find rotation of delta x delta square projected on unit sphere
delta = 0.1;

% radius of sphere
r = 1;

% construct square
points = [-delta -delta;
    -delta delta;
    delta delta;
    delta -delta];
points_old = points;

% calculate 2d rotation matrix
dtheta = dtheta * (pi/180);
R_2d = [cos(dtheta) -sin(dtheta);
    sin(dtheta) cos(dtheta)];

% rotate square in 3d
points = (R_2d*points')';

% translate points in 2d
points(:,1) = points(:,1) + (dxpixels/npixels)*(2*xmax);
points(:,2) = points(:,2) + (dypixels/npixels)*(2*ymax);

% project points onto unit sphere
points = [points sqrt(r^2-points(:,1).^2-points(:,2).^2)];
points_old = [points_old sqrt(r^2-points_old(:,1).^2-points_old(:,2).^2)];

% calculate covariance matrix
C = points_old'*points;

% do svd
[u, s, v] = svd(C);

d = sign(det(C));
I = eye(3);
I(end,end) = d;

R = v * I * u';


    


