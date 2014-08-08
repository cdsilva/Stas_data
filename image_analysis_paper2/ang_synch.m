function R_opt = ang_synch(R, dim)

m = size(R,1) / dim;
if mod(m, 1) ~= 0
    disp('dim is not compatable with size of R')
    return
end


[V, D] = eigs(R, dim);
[~, ind] = sort(abs(diag(D)), 'descend');
V = V(:, ind);
D = D(ind, ind);
for i=1:dim
    V(:,i) = V(:,i) / norm(V(:,i));
end
if det(V(1:dim,1:dim)*D(1:dim,1:dim)*V(1:dim,1:dim)') < 0
    V(:,dim) = -V(:,dim);
end

R_opt = V;
for i=1:m
    [u, s, v] = svd(R_opt(dim*(i-1)+1:dim*i,:));
    R_opt(dim*(i-1)+1:dim*i,:) = v * u';
end
for i=m:-1:1
    R_opt(dim*(i-1)+1:dim*i,:) = R_opt(1:dim,:)' * R_opt(dim*(i-1)+1:dim*i,:);
end


% R_opt = V*D*V(1:dim,:)';
% 
% for i=1:m
%    R_opt(dim*(i-1)+1:dim*i,:) = R_opt(dim*(i-1)+1:dim*i,:) /(det(R_opt(dim*(i-1)+1:dim*i,:))^(1/dim));
% end

% for i=1:m
%     det(R_opt(dim*(i-1)+1:dim*i,:))
%     R_opt(dim*(i-1)+1:dim*i,:)'*R_opt(dim*(i-1)+1:dim*i,:)
% end

