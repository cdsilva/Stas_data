function [I, V, D] = unscramble_dmaps(data, eps)

W = squareform(pdist(data));
W = W.^2;

%% run dmaps

[V, D]=dmaps(W,eps,6);

%% find time embedding from first (non-trivial) dmaps coordinate
[~,I] = sort(V(:,2));

if sum(data(I(1),:)) > sum(data(I(end),:))
    I = flipud(I);
end