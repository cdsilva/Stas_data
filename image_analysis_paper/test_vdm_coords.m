clear all
close all

dim = 3;
n = 100;

R_true = zeros(n*dim, dim);
R_true(1:dim,:) = eye(dim);
for i=dim+1:dim:dim*(n-1)+1
    [u, s, v] = svd(rand(dim));
    if det(u*v') < 0
        u(:,end) = -u(:,end);
    end
    R_true(i:i+dim-1,:) = u * v';
end
R = R_true * R_true';

W = inf(n);
W(1,1) = 0;
W(1,2) = 1;
for i=2:n-1
    W(i, i) = 0;
    W(i, i+1) = 1;
    W(i, i-1) = 1;
end
W(n, n) = 0;
W(n, n-1) = 1;

eps = 1;
neigs = 10;

[R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs);
    
    
figure;
plot(embed_coord(:, 7))
    