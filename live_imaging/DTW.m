function [D, warp_path] = DTW(d)
%% d is the pairwise distance matrix
%% D is the time warp matrix (cumulative distance over paths)
%% warp_path contains the optimal time warp path (warp_path(i,j) = 1 if points i and j are connected/linked)

[m, n] = size(d);

D = zeros(size(d));

D(1,:) = cumsum(d(1,:));

for i=2:m
    j = 1;
    D(i,j) = d(i,j) + D(i-1, j);
    for j=2:n
        D(i,j) = d(i,j) + min([D(i-1, j), D(i, j-1), D(i-1, j-1)]);
    end
end

%%
warp_path = zeros(size(d));

i = m;
j = n;
warp_path(i,j) = 1;

while i > 1 || j > 1
    if i > 1
        i_tmp = i-1;
        j_tmp = j;
        if j > 1
            if D(i-1,j-1) < D(i_tmp, j_tmp)
                i_tmp = i-1;
                j_tmp = j-1;
            end
            if D(i,j-1) < D(i_tmp, j_tmp)
                i_tmp = i;
                j_tmp = j-1;
            end
        end
    else
        i_tmp = i;
        j_tmp = j-1;
    end
    i = i_tmp;
    j = j_tmp;
    warp_path(i,j) = 1;
end
