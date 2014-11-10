function [cp1, cp_deltas] = optical_flow_warping(im1, im2)

npixels = size(im1, 1);
filt_size = round(0.02*npixels);
[X, Y] = meshgrid(1:npixels,1:npixels);

% im1 = medfilt2(adapthisteq(im1), [filt_size filt_size]);
im1 = medfilt2(im1, [filt_size filt_size]);
im2 = medfilt2(im2, [filt_size filt_size]);

%%
[u, v] = HS(im1, im2);

%%

navg = 20;
stride = npixels/navg;

avg_X = zeros(navg);
avg_Y = zeros(navg);
avg_u = zeros(navg);
avg_v = zeros(navg);

for i=1:navg
    for j=1:navg
        avg_X(i,j) = mean(mean(X((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_Y(i,j) = mean(mean(Y((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_u(i,j) = mean(mean(u((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
        avg_v(i,j) = mean(mean(v((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride)));
    end
end
    
%%
cp1 = [reshape(avg_X, [], 1) reshape(avg_Y, [], 1)];
cp_deltas = [reshape(avg_u, [], 1) reshape(avg_v, [], 1)];

% %%
% 
% sigma = 0.1;
% im1_warped = warp_image_rbf(im1, [reshape(avg_X, [], 1) reshape(avg_Y, [], 1)], [reshape(avg_X+avg_u, [], 1) reshape(avg_Y+avg_v, [], 1)], sigma);

