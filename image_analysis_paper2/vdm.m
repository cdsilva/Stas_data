% function [R_opt, embed_coord, embed_idx, D] = vdm(R, W, eps, neigs)
function [R_opt, embed_coord, D2, D] = vdm(R, W, eps, neigs)

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
    %     R_opt(dim*(i-1)+1:dim*i,:) = v * u';
    R_opt(dim*(i-1)+1:dim*i,:) = u * v';
end
for i=n:-1:1
    R_opt(dim*(i-1)+1:dim*i,:) = R_opt(1:dim,:)' * R_opt(dim*(i-1)+1:dim*i,:);
end

% embed_coord = zeros(n, neigs*(neigs-1)/2+neigs);
% embed_idx = zeros(2, neigs*(neigs-1)/2+neigs);
% curr_idx = 1;
% for i=1:neigs
%     for j=1:i
%         embed_coord(:, curr_idx) = sum(reshape(V(:,i),dim, []).*reshape(V(:,j),dim, []))';
%         embed_idx(1, curr_idx) = i;
%         embed_idx(2, curr_idx) = j;
%         curr_idx = curr_idx + 1;
%     end
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




