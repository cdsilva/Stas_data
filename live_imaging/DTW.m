function D = DTW(d)

[m, n] = size(d);

D = zeros(size(d));

D(1, 2:end) = inf;
D(2:end, 1) = inf;

for i=2:n
    for j=2:(m-i)
        D(i,j) = d(i,j) + min([D(i-1, j) D(i-1,j-1) D(i,j-1)]);
    end
end


