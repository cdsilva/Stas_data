clear all
close all

n = 20;
theta = linspace(0, 3*pi/2, n);

x = cos(theta);
y = sin(theta);

%x = -1:0.1:1;
%y = x.^3;

cm = gray;
cm = flipud(cm);


W = squareform(pdist([x' y'])).^2;
eps = median(W(:))/2;

A = exp(-W/eps);
%A = diag(1./sum(A))*A;

figure;
plot(x, y, '.', 'markersize', 10)


figure;
for i=1:n
    for j=1:n-1
        idx = ceil(size(cm, 1) / max(A(:)) * A(i,j));
        plot([x(i) x(j)], [y(i) y(j)], 'color', cm(idx, :), 'linewidth', 1);
        hold on
    end
end     
plot(x, y, '.', 'markersize', 20)

[V, D] = dmaps(W, eps, 2);
figure;
scatter(x, y, 500, V(:,2), '.')

