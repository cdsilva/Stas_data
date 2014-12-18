function r = compute_ranks(x)
% r = compute_ranks(x)
% compute the rank of each of the entries in x
% ranks are returned in r

r = zeros(size(x));

m = length(x);

[~, ind] = sort(x);

r(ind) = 1:m;
