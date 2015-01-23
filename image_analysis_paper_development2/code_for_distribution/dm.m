function [embed_coord, D2] = dm(W, eps_scale, ncomps)

eps = median(W(:)) * eps_scale;

neigs = ncomps + 1;

% calculate kernel of distances
W2 = exp(-(W/eps).^2);
W2 = diag(1./sum(W2)) * W2;

% compute eigenvectors
[V, D] = eigs(W2, neigs);

[~, I] = sort(abs(diag(D)), 'descend');

V = V(:,I);
D = D(I,I);

embed_coord = V(:, 2:end);
D2 = diag(D(2:end, 2:end));


