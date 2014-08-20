clear all
close all

%%


nimages_all = 132;
image_dir = '../membrane_pictures/14_0501_dpERK_late';
image_name = 'emb';

shift_max = 0.1;

neigs = 10;
alpha = 0;

channel = 1;
nchannels = 3;

image_set_raw_all = zeros(1024, 1024, nchannels, nimages_all, 'uint8');

%%

for i=1:nimages_all
    im_tmp = imread(sprintf('%s/%s%02d.tif', image_dir, image_name, i));
    image_set_raw_all(:, :, :, i) = im_tmp;
end

%%

% npixels_vec = [50 100 200];
npixels_vec = [50 100 200];
nimages_vec = [20 50 100 120];
nrot_vec = [10 20 30 40];
% nshifts_vec = [0 5 10];
nshifts_vec = [0];

timing = zeros(length(npixels_vec), length(npixels_vec), length(nrot_vec), length(nshifts_vec));

%%

for i=1:length(npixels_vec)
    image_set_all = imresize(image_set_raw_all, [npixels_vec(i) npixels_vec(i)]);
    for j=1:length(nimages_vec)
        nimages = nimages_vec(j);
        image_set = image_set_all(:, :, :, 1:nimages);
        for k=1:length(nrot_vec)
            nrot = nrot_vec(k);
            for l=1:length(nshifts_vec)
                nshifts = nshifts_vec(l);
                
                tic
                
                %
                [R, W] = compute_pairwise_alignments(image_set, nrot, nshifts, shift_max);
                dim = size(R, 1) / nimages;
                
                %
                W2 = W.^2;
                eps = median(W2(:))/10;
                [R_opt, embed_coord, D2, D] = vdm(R, W2, eps, neigs, alpha);
                
                timing(i,j,k,l) = toc;
            end
        end
    end
end

%%

figure;
for i=1:length(npixels_vec)
    for j=1:length(nimages_vec)
        % for k=1:length(nrot_vec)
        for l=1:length(nshifts_vec)
            loglog(nrot_vec, squeeze(timing(i,j,:,l)), '-o')
            hold on
        end
    end
end
loglog([10 50], [10 50], '-r')
xlabel('number of rotations')
ylabel('time (sec)')

figure;
for i=1:length(npixels_vec)
    % for j=1:length(nimages_vec)
    for k=1:length(nrot_vec)
        for l=1:length(nshifts_vec)
            loglog(nimages_vec, squeeze(timing(i,:,k,l)),'-o')
            hold on
        end
    end
end
loglog([10 100], ([10 100].^2)/100, '-r')

xlabel('number of images')
ylabel('time (sec)')

figure;
% for i=1:length(npixels_vec)
    for j=1:length(nimages_vec)
    for k=1:length(nrot_vec)
        for l=1:length(nshifts_vec)
            loglog(npixels_vec, squeeze(timing(:,j,k,l)),'-o')
            hold on
        end
    end
end
loglog([10 100], ([10 100].^2)/100, '-r')
xlabel('number of images')
ylabel('time (sec)')
