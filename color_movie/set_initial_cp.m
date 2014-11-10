function cp = set_initial_cp(im1, ncp)

npixels = size(im1, 1);

im1 = adapthisteq(im1);

thres = 50;
idx = find(im1 > thres);
cpi = zeros(ncp, 2);
[cpi(:,1), cpi(:,2)] = ind2sub(npixels, randsample(idx, ncp));

% cpi = randi(npixels, [ncp, 2]);

options = optimset('MaxFunEvals',1e6);
cp = fminsearch(@(x) total_cost(x, im1), cpi, options);

function y = total_cost(cp, im1)

repel_force = 1e-1;
sigma = 1e-2;

npixels = size(im1, 1);

y = repel_force * sum(exp(-pdist(cp).^2/(sigma^2*npixels^2))) - sum(interp2(1:npixels,1:npixels,double(im1),cp(:,1),cp(:,2)));







