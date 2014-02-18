function [R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs)

dim = size(R,1) / size(W,1);

[n, ~] = size(W);

if mod(dim, 1) ~= 0
    disp('sizes of R and W are not compatible')
    return
end

R2 = zeros(size(R));
W2 = exp(-W/eps);
W2 = diag(1./sum(W2)) * W2;

for i=1:n
    for j=1:n
        R2(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j) = W2(i,j) * R(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j);
    end
end

[V, D] = eigs(R2, neigs);
[~, ind] = sort(abs(diag(D)), 'descend');
V = V(:, ind);
D = D(ind, ind);
for i=1:neigs
    V(:,i) = V(:,i) / norm(V(:,i));
end
if det(V(1:dim,1:dim)*D(1:dim,1:dim)*V(1:dim,1:dim)') < 0
    V(:,dim) = -V(:,dim);
end

R_opt = V(:,1:dim);
for i=1:n
    [u, s, v] = svd(R_opt(dim*(i-1)+1:dim*i,:));
    R_opt(dim*(i-1)+1:dim*i,:) = u * v';
end
R_opt = R_opt * R_opt(1:dim,:)';

% 
% R_opt = V(:,1:dim)*D(1:dim,1:dim)*V(1:dim,1:dim)';
% 
% for i=1:n
%     R_opt(dim*(i-1)+1:dim*i,:) = R_opt(dim*(i-1)+1:dim*i,:) /(det(R_opt(dim*(i-1)+1:dim*i,:))^(1/dim));
% end

embed_coord = zeros(n, neigs*(neigs-1)/2+neigs);
embed_idx = zeros(2, neigs*(neigs-1)/2+neigs);
curr_idx = 1;
for i=1:neigs
    for j=1:i
        embed_coord(:, curr_idx) = sum(reshape(V(:,i),dim, []).*reshape(V(:,j),dim, []))';
        embed_idx(1, curr_idx) = i;
        embed_idx(2, curr_idx) = j;
        curr_idx = curr_idx + 1;
    end
end
        