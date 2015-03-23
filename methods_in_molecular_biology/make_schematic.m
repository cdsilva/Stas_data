clear all
close all


%%

% n1 = 100;
% n2 = 100;
% 
% x = [linspace(0, 1, n1)'; cos(linspace(0, pi, n2))'];
% y = [zeros(n1, 1); sin(linspace(0, pi, n2))'];
theta = linspace(0, 2*pi, 200)';

x = theta.*cos(theta);
y = theta.*sin(theta);

t = (1:length(x))';

figure;
scatter(x,y,50, t, '.')
colorbar
axis equal

%%

data = [x y];
mdl = fitlm(data,t,'linear');

[pred_time, pred_time_int] = predict(mdl, data);

figure;
plot(t, pred_time, '.')

figure;
scatter(x,y,50, pred_time, '.')
colorbar
axis equal

[X, Y] = meshgrid(linspace(min(x),max(x),100), linspace(min(y),max(y),100));
[pred_time, pred_time_int] = predict(mdl, [X(:) Y(:)]);

figure;
imagesc(reshape(pred_time, size(X)))
colorbar