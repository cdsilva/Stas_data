function r = compute_ranks(x)

r = zeros(size(x));

m = length(x);

[~, ind] = sort(x);

r(ind) = 1:m;
