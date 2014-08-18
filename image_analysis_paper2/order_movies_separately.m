clear all
close all

%%

npixels = 100;

dt = 0.5;

nmovies = 6;
make_subplot = @(k) subplot(2, 3, k);

% idx_start = [21 31 27 24 18 17];
idx_start = ones(1, 6);
idx_end = [56 60 47 44 46 42];

f1 = figure;
f2 = figure;

for k=1:nmovies
    load(sprintf('movie%d.mat', k));
%     idx = idx_start(k):idx_end(k);
    idx = idx_end(k)-41:idx_end(k);
    
    time = time * dt;
    
    nimages = length(time);
    
    image_set = zeros(npixels, npixels, nimages, 'uint8');
    channel = 1;
    
    filt = fspecial('disk');
    
    for i=1:nimages
        image = images(:, :, channel, i);
        image = imresize(image, [npixels npixels]);
        %     image = adapthisteq(image, 'cliplimit',0.005);
        %     image = adapthisteq(image, 'numtiles', [16 16]);
        image = crop_image(image, channel);
        image = padarray(image, [10 10]);
        image = imresize(image, [npixels npixels]);
        %     image = adapthisteq(image, 'cliplimit',0.05);
        image = adapthisteq(image);
        %     image = medfilt2(image);
        image = imfilter(image, filt, 'replicate');
        
        image_set(:, :, i) = image;
    end
    
    %%
    image_set = image_set(:, :, idx);
    time = time(idx);
    nimages = length(time);
    
    %%
    
    rng(21*k);
    theta = 360*rand(nimages, 1);
    idx = randperm(nimages);
    
    for i=1:nimages
        image_set(:, :, i) = imrotate(image_set(:, :, i), theta(i), 'crop');
    end
    
    image_set = image_set(:, :, idx);
    time = time(idx);
    theta = theta(idx);
    
    %%
    
    nrot = 40;
    nshifts = 0;
    shift_max = 0.1;
    
    [R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
    
    dim = size(R, 1) / nimages;
    
    %%
    W2 = W.^2;
    eps = median(W2(:))/10;
    neigs = 10;
    alpha = 0;
    [R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);
    
    %%
    theta_opt = zeros(size(theta));
    for i=1:nimages
        R_tmp = R_opt(dim*(i-1)+1:dim*i, :);
        theta_opt(i) = atan2d(R_tmp(2,1), R_tmp(1,1));
    end
    
    %%
    figure(f1);
    make_subplot(k)
    plot(mod(theta, 360), mod(theta_opt, 360), '.')
    xlabel('true angle')
    ylabel('recovered angle')
    title(sprintf('average error = %2.2f', std(mod(theta - theta_opt, 360))))
    
    %%
    if corr(time, embed_coord(:,1)) < 0
        embed_coord(:,1) = -embed_coord(:,1);
    end
    
    
    figure(f2);
    make_subplot(k)
    plot(compute_ranks(time), compute_ranks(embed_coord(:,1)),'.')
    xlabel('rank from time')
    ylabel('rank from vdm')
    title(sprintf('rank corr = %2.4f', corr(compute_ranks(time), compute_ranks(embed_coord(:,1)))))
    
    
end