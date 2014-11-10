function [cp1, cp2] = match_contours(B1, B2)

D = pdist2(B1{1}, B2{1});
[~, I] = min(D, [], 2);

cp1 = B1{1};
cp2 = B2{1}(I,:);

%%
D = pdist2(B1{2}, B2{2});
[~, I] = min(D, [], 2);

cp1 = [cp1; B1{2}];
cp2 = [cp2; B2{2}(I,:)];

cp1 = fliplr(cp1);
cp2 = fliplr(cp2);
