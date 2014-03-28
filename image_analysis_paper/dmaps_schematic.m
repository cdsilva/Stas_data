clear all
close all

res = '-r300';
fmt = '-djpeg';

fontsize = 10;

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
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
for i=1:n
    for j=1:n-1
        idx = ceil(size(cm, 1) / max(A(:)) * A(i,j));
        plot([x(i) x(j)], [y(i) y(j)], 'color', cm(idx, :), 'linewidth', 0.5);
        hold on
    end
end     
scatter(x, y, 200, 'k', '.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
box off
xlabel('x', 'fontsize', fontsize)
ylabel('y', 'fontsize', fontsize)
%print(sprintf('%s/dmaps_schematic_edges', im_save_dir), fmt, res);
saveas(gcf, sprintf('%s/fig4a', im_save_dir), 'epsc');


[V, D] = dmaps(W, eps, 2);
figure;
set(gcf, 'paperunits', 'centimeters')
set(gcf, 'papersize', [4 4])
set(gcf, 'paperposition',[0 0 4 4])
scatter(x, y, 200, V(:,2), '.')
cm_red = gray;
cm_red(:,2) = 0;
cm_red(:,3) = 0;
colormap(cm_red)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('x', 'fontsize', fontsize)
ylabel('y', 'fontsize', fontsize)
%print(sprintf('%s/dmaps_schematic_color', im_save_dir), fmt, res);
saveas(gcf, sprintf('%s/fig4b', im_save_dir), 'epsc');
