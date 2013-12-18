function [R, W, theta] = align_data(data)

[m, n] = size(data);

R = zeros(2*m);
theta = zeros(m);
W = zeros(m);

for i=1:m
    all_rot = gallery('circul',data(i,:));
    for j=1:i-1
        [min_dist, ind] = min(sum((all_rot - ones(n,1)*data(j,:)).^2, 2));
        ind = ind - 1;
        theta(i,j) = 2 * pi * ind / n;
        theta(j,i) = -theta(i,j);
        
        R(2*i-1:2*i,2*j-1:2*j) = [cos(theta(i,j)) -sin(theta(i,j)); sin(theta(i,j)) cos(theta(i,j))];
        R(2*j-1:2*j,2*i-1:2*i) = R(2*i-1:2*i,2*j-1:2*j)';
        
        W(i,j) = min_dist;
        W(j,i) = min_dist;
    end
end


