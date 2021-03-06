% function [R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs)
function [R_opt, embed_coord, D2, D] = vdm(R, W, eps, neigs, alpha)

dim = size(R,1) / size(W,1);

[n, ~] = size(W);

if mod(dim, 1) ~= 0
    disp('sizes of R and W are not compatible')
    return
end

R2 = zeros(size(R));
W2 = exp(-W/eps);
for i=1:n
    W2(i, i) = 0;
end

D = sum(W2);

if alpha > 0
    W2 = diag(D.^(-alpha)) * W2 * diag(D.^(-alpha));
end
W2 = diag(1./sum(W2)) * W2;

for i=1:n
    for j=1:n
        R2(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j) = W2(i,j) * R(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j);
    end
end

[V, D] = eigs(R2, neigs);
imag_tol = 1e-4;
count = 0;
count_max = 100;
while norm(imag(D)) < imag_tol && norm(imag(V)) > imag_tol && count < count_max
    [V, D] = eigs(R2, neigs);
    count = count + 1;
end    
if count == count_max || norm(imag(D)) > imag_tol
    disp('ERROR: imaginary eigenvectors')
end

[~, ind] = sort(abs(diag(D)), 'descend');
V = V(:, ind);
D = D(ind, ind);
for i=1:neigs
    V(:,i) = V(:,i) / norm(V(:,i));
end
if det(V(1:dim,1:dim)) < 0
    V(:,dim) = -V(:,dim);
end

% R_opt = V(:,1:dim) * V(1:dim,1:dim)';
R_opt = V(:,1:dim);
for i=1:n
    [u, s, v] = svd(R_opt(dim*(i-1)+1:dim*i,:));
    %     R_opt(dim*(i-1)+1:dim*i,:) = v * u';
    R_opt(dim*(i-1)+1:dim*i,:) = u * v';
end
R_opt = R_opt * R_opt(1:dim,1:dim)';
% for i=n:-1:1
%     R_opt(dim*(i-1)+1:dim*i,:) = R_opt(1:dim,:)' * R_opt(dim*(i-1)+1:dim*i,:);
% end

embed_coord = zeros(n, floor(neigs/dim)-1);
D2 = zeros(floor(neigs/dim)-1, 1);
for i=1:size(embed_coord, 2);
    var_coord = 0;
    for j1=1:dim
        for j2=dim*i+1:dim*(i+1)
            embed_coord_tmp = sum(reshape(V(:,j1),dim, []).*reshape(V(:,j2),dim, []))';
            var_tmp = var(embed_coord_tmp);
            if var_tmp > var_coord
                var_coord = var_tmp;
                embed_coord(:, i) = embed_coord_tmp;
                D2(i) = D(j1,j1) * D(j2,j2);
            end
        end
    end
end




