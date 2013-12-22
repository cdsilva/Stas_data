function [I, coeff, V_PCA, D_PCA] = unscramble_pca(data, symmetrize)
% computes ordering of data using PCA
% data is the data matrix (each row is a data point)
% symmetrize is a boolean argument; true means to first symmetrize the data
% I is the indices that sort the data by the projection coeff. onto the
% first PC
% coeff are the projection coefficients for the data onto the PCs (each row
% is a data point, each column is a PC)
% V_PCA stores the PCA (one PC per column)
% D_PCA is a diaganol matrix with the eigenvalues for each PC

[m, n] = size(data);

if nargin == 1 || symmetrize == false
    [V_PCA, D_PCA] = PCA(data, n);
else
    data2 = [data; data(:,end:-1:1)];
    [V_PCA, D_PCA] = PCA(data2, n);
end

% mean - center
mean_data = mean(data, 1);
data = data - ones(m, 1) * mean_data;

coeff = data * V_PCA;
[~, I] = sort(coeff(:,1));

if sum(data(I(1),:)) > sum(data(I(end),:))
    I = flipud(I);
    V_PCA = -V_PCA;
    coeff = -coeff;
end
