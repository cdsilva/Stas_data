clear all
close all

%%
% <latex>
% First, we create an artifical data set of $n=7$ data points which lie 
% on a one-dimensional curve in two dimensions
% </latex>

n = 7;

theta = linspace(0, 3*pi/2, n);

x = [0.9*cos(theta); 0.9*sin(theta)]

figure;
scatter(x(1,:),x(2,:),50, 'k', '.');
xlabel('x(1)');
ylabel('x(2)');

%%
% <latex>
% Then we compute the $n \times n$ distance matrix between data points
% </latex>

dist_matrix = dist(x)

%%
% <latex>
% We choose $\epsilon$, the kernel scale, as half the median of the
% pairwise distances
% </latex>

eps = median(dist_matrix(:))/2

%%
% <latex>
% We then compute the $n \times n$ matrix $W$, with 
% \[ W_{ij} = \exp \left(-\frac{ \| x_i - x_j \|^2}{\epsilon^2} \right)\]
% </latex>

W = exp(-dist_matrix.^2/eps.^2)


%%
% <latex>
% We then compute the $n \times n$ diaganol matrix $D$, with 
% \[ D_{ii} = \sum_{j=1}^n W_{ij} \]
% </latex>

D = diag(sum(W))

%%
% <latex>
% We then compute the $n \times n$ matrix $A = D^{-1}W$
% </latex>

A = inv(D)*W

%%
% <latex>
% We then compute the top eigenvectors and eigenvalues of $A$, and order
% them by the magnitude of the eigenvalues
% </latex>

[evecs, evals] = eig(A);
[~, I] = sort(abs(diag(evals)), 'descend');

evecs = evecs(:,I)
evals = evals(I,I)

%%
% <latex>
% The first eigenvector is a trivial constant eigenvector with eigenvalue
% 1. The second eigenvector provides an ordering of the data. If we color
% our data by the value of this second eigenvector, we see it is one-to-one
% with the arclength along the curve. Furthermore, if we sorted the data by
% the corresponding entries of the second eigenvector, we would order the
% data along this one-dimensional nonlinear curve. 
% </latex>

figure;
scatter(x(1,:),x(2,:),50, evecs(:,2), '.');
xlabel('x(1)');
ylabel('x(2)');



