clear all
close all

%%

npixels = 100;

dt = 0.5;

idx_end = [56 60 47 44 46 42];
idx_start = idx_end - 41;

k = 2;
rng(123);

load(sprintf('movie%d.mat', k));
idx = idx_start(k):idx_end(k);

time_all = time * dt;

nimages_all = length(time_all);

image_set_all = zeros(npixels, npixels, nimages_all, 'uint8');
channel = 1;

filt = fspecial('disk');

for i=1:nimages_all
    image = images(:, :, channel, i);
    image = imresize(image, [npixels npixels]);
    
    image = crop_image(image, channel);
    image = padarray(image, [10 10]);
    image = imresize(image, [npixels npixels]);
    
    image = adapthisteq(image);
    
    image = imfilter(image, filt, 'replicate');
    
    image_set_all(:, :, i) = image;
end

%%
image_set_all = image_set_all(:, :, idx);
time_all = time_all(idx);
nimages_all = length(time_all);

%%

nreps = 25;

sample_vec = round(linspace(10, nimages_all, 10));
theta_err = zeros(length(sample_vec), nreps);
rank_corr = zeros(length(sample_vec), nreps);

for j = 1:length(sample_vec)
    nimages = sample_vec(j);

    j
    
    for j2 = 1:nreps
        
        sample_idx = randi(nimages_all, nimages, 1);
        
        image_set = image_set_all(:, :, sample_idx);
        time = time_all(sample_idx);
        
        
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
        
        theta_err(j, j2) = std(mod(theta-theta_opt, 360));
        
        if theta_err(j, j2) > 20
            theta_err(j, j2) = std(mod(theta-theta_opt+180, 360)-180);
        end
        
        %%
        if corr(time, embed_coord(:,1), 'type','spearman') < 0
            embed_coord(:,1) = -embed_coord(:,1);
        end
        
        rank_corr(j, j2) = corr(compute_ranks(time), compute_ranks(embed_coord(:,1)));
    end
end

%%
figure;
plot(sample_vec, mean(theta_err, 2), '.')

figure;
plot(sample_vec, mean(rank_corr, 2), '.')
xlabel('number of images')
ylabel('average rank correlation coefficient')

figure;
errorbar(sample_vec, mean(theta_err, 2), std(theta_err,[], 2))

figure;
errorbar(sample_vec, mean(rank_corr, 2), std(rank_corr,[], 2))
