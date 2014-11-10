% clear all
close all

%% read in images

% npixels = 500;
% 
% % gastrulation
% [images1, time1] = load_fixed_images1(npixels);
% time1 = time1 + 50;
% nimages1 = length(time1);
% 
% % cellularization
% [images2, time2] = load_fixed_images2(npixels);
% nimages2 = length(time2);
% 
% % ind cell + gast
% [images3, time3] = load_fixed_images3(npixels);
% nimages3 = length(time3);
% 
% % movie
% [images_movie, time_movie] = load_movie(npixels);
% nimages_movie = length(time_movie);

%%
nchannels = 4;
nuclei_channel = 1;
dpERK_channel = 2;
Dl_channel = 3;
ind_channel = 4;

all_nimages = nimages1 + nimages2 + nimages3;

all_images = zeros(npixels, npixels, nchannels, all_nimages, 'uint8');
has_channels = zeros(all_nimages, nchannels);

% images1: nuclei (red), dpERK (green), and Dl (blue)
all_images(:, :, nuclei_channel, 1:nimages1) = images1(:,:,1,:);
all_images(:, :, dpERK_channel, 1:nimages1) = images1(:,:,2,:);
all_images(:, :, Dl_channel, 1:nimages1) = images1(:,:,3,:);
has_channels(1:nimages1, [nuclei_channel dpERK_channel Dl_channel]) = 1;

% images2: nuclei (red), dpERK (green), and Dl (blue)
all_images(:, :, nuclei_channel, nimages1+(1:nimages2)) = images2(:,:,1,:);
all_images(:, :, dpERK_channel, nimages1+(1:nimages2)) = images2(:,:,2,:);
all_images(:, :, Dl_channel, nimages1+(1:nimages2)) = images2(:,:,3,:);
has_channels(nimages1+(1:nimages2), [nuclei_channel dpERK_channel Dl_channel]) = 1;

% images3: dpERK (red), ind (green), and nuclei (blue)
all_images(:, :, nuclei_channel, nimages1+nimages2+(1:nimages3)) = images3(:,:,3,:);
all_images(:, :, dpERK_channel, nimages1+nimages2+(1:nimages3)) = images3(:,:,1,:);
all_images(:, :, ind_channel, nimages1+nimages2+(1:nimages3)) = images3(:,:,2,:);
has_channels(nimages1+nimages2+(1:nimages3), [nuclei_channel dpERK_channel ind_channel]) = 1;

all_time = [time1; time2; time3];

%%

images_movie_colored = zeros(npixels, npixels, nchannels, nimages_movie, 'uint8');
kernel_scale = 1;
weight_thres = 1e-2;
for i=1:nimages_movie
% for i=[120]
    i
    weights = exp(-(all_time - time_movie(i)).^2/kernel_scale.^2);
    weights = weights / sum(weights);
    
    total_weights = zeros(nchannels, 1);
    
    for j=find(weights > weight_thres)'
        
%         [cp1, cp_deltas] = optical_flow_warping(all_images(:,:,nuclei_channel, j), images_movie(:,:,nuclei_channel,i));
%         im1_warped = warp_image_rbf(all_images(:,:,:,j), cp1, cp1+cp_deltas, 0.05);
        [cp1, cp_deltas] = optical_flow_warping(images_movie(:,:,nuclei_channel,i),adapthisteq(all_images(:,:,nuclei_channel, j)));
        im1_warped = warp_image_rbf2(all_images(:,:,:,j), cp1, cp_deltas, 0.05);
        
%         figure;
%         imshow(imabsdiff(images_movie(:,:,nuclei_channel,i),im1_warped(:,:,nuclei_channel)))
%         subplot(2,2,1)
%         imshow(images_movie(:,:,nuclei_channel,i))
%         subplot(2,2,2)
%         imshow(all_images(:,:,nuclei_channel, j))
%         subplot(2,2,3)
%         imshow(images_movie(:,:,nuclei_channel,i))
%         subplot(2,2,4)
%         imshow(im1_warped(:,:,nuclei_channel))
        
        for k=2:nchannels
            if has_channels(j, k) == 1
                images_movie_colored(:,:,k,i) = imlincomb(1, images_movie_colored(:,:,k,i), weights(j), im1_warped(:,:,k));
                total_weights(k) = total_weights(k) + weights(j);
            end
        end
    end
    for k=2:nchannels
        images_movie_colored(:,:,k,i) = imdivide(images_movie_colored(:,:,k,i), total_weights(k));
    end
    images_movie_colored(:,:,nuclei_channel,i) = images_movie(:,:,nuclei_channel,i);
    
end

%%

combine_channels = @(x) imlincomb(0.25,repmat(x(:,:,nuclei_channel), [1 1 3]), 1, x(:,:,2:end));

writerObj = VideoWriter('colored_movie_warped.avi');
writerObj.FrameRate = 10;
open(writerObj);

figure;

imshow(combine_channels(images_movie_colored(:,:,:,1)), 'initialmagnification','fit','border','tight')
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
clf

for i=1:nimages_movie
    imshow(combine_channels(images_movie_colored(:,:,:,i)))
    frame = getframe;
    writeVideo(writerObj,frame);
    clf
end

close(writerObj);


