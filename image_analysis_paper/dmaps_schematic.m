clear all
close all

res = '-r300';
fmt = '-djpeg';

im_save_dir = 'paper_figures';
%%
n = 15;
theta = linspace(0, 3*pi/2, n);

x = 0.9*cos(theta);
y = 0.9*sin(theta);

%x = -1:0.1:1;
%y = x.^3;

cm = gray;
cm = flipud(cm);


W = squareform(pdist([x' y'])).^2;
eps = median(W(:))/2;

A = exp(-W/eps);
%A = diag(1./sum(A))*A;


%%

figure;
for i=1:n
    for j=1:n-1
        idx = ceil(size(cm, 1) / max(A(:)) * A(i,j));
        plot([x(i) x(j)], [y(i) y(j)], 'color', cm(idx, :), 'linewidth', 2);
        hold on
    end
end     
plot(x, y, '.', 'markersize', 50)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
%axis off
print(sprintf('%s/dmaps_schematic_edges', im_save_dir), fmt, res);


[V, D] = dmaps(W, eps, 2);
figure;
scatter(x, y, 2000, V(:,2), '.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print(sprintf('%s/dmaps_schematic_color', im_save_dir), fmt, res);

