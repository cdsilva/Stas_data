function [R, W] = compute_pairwise_alignments_1d_color(dpERK)

dim = 2;

[m, n, p] = size(dpERK);

R = zeros(dim*m);
W = inf(m);

thetas = zeros(m);

for i=1:m
    imagei = dpERK(i, :, :);
    for k=1:n
        imagei_tmp = circshift(imagei, [0 k 0]);
        for j=1:i-1
            imagej = dpERK(j, :, :);
            d2 = norm(imagej(:) - imagei_tmp(:))^2;
            if d2 < W(i, j)
                W(i,j) = d2;
                thetas(j, i) = k/n * 2 * pi;
            end
        end
    end
end


for i=1:m
    for j=1:i-1
        W(j,i) = W(i,j);
        R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i) = calc_rot_matrix(thetas(j, i));
        R(dim*(i-1)+1:dim*i, dim*(j-1)+1:dim*j) = R(dim*(j-1)+1:dim*j, dim*(i-1)+1:dim*i)';
    end
end

function R = calc_rot_matrix(theta)

R = [cos(theta) sin(theta);
    -sin(theta) cos(theta)];





