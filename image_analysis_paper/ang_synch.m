function R_opt = ang_synch(R, dim)

m = size(R,1) / dim;
if mod(m, 1) ~= 0
    disp('dim is not compatable with size of R')
    return
end


[V, D] = eigs(R, dim);
if det(V(1:dim,:)*D*V(1:dim,:)') < 0
    V(:,end) = -V(:,end);
end

R_opt = V*D*V(1:dim,:)';

for i=1:m
    R_opt(dim*(i-1)+1:dim*i,:) = R_opt(dim*(i-1)+1:dim*i,:) /(det(R_opt(dim*(i-1)+1:dim*i,:))^(1/dim));
end