function theta = R_to_theta(R_opt)
% theta = R_to_theta(R_opt)
% converts the rotation matrices to angles
% R_opt is a 2nx2 matrix, where each 2x2 block corresponds to a rotation
% matrix
% theta is a n-long vector containing the angles correponding to each
% rotation matrix in R_opt

[m, dim] = size(R_opt);

if dim ~= 2
    disp('ERROR: R_opt is not of the correct dimensions');
    return
end

n = m / dim;
if mod(n, 1) ~= 0
    disp('ERROR: R_opt is not of the correct dimensions');
    return
end

theta = zeros(n, 1);
for i=1:n
    R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
    theta(i) = atan2d(R_tmp(2,1), R_tmp(1,1));
end


