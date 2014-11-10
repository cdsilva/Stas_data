clear all
close all

%% read in images

npixels = 500;

% gastrulation
[images1, time1] = load_fixed_images1(npixels);
time1 = time1 + 50;
nimages1 = length(time1);
nuclei_contours1 = cell(nimages1, 1);
for i=1:nimages1
    B = find_active_contours(images1(:,:,1,i));
    nuclei_contours1{i} = B;
end

% cellularization
[images2, time2] = load_fixed_images2(npixels);
nimages2 = length(time2);
nuclei_contours2 = cell(nimages2, 1);
for i=1:nimages2
    B = find_active_contours(images2(:,:,1,i));
    nuclei_contours2{i} = B;
end

% ind cell + gast
[images3, time3] = load_fixed_images3(npixels);
nimages3 = length(time3);
nuclei_contours3 = cell(nimages3, 1);
for i=1:nimages3
    B = find_active_contours(images3(:,:,3,i));
    nuclei_contours3{i} = B;
end

% movie
[images_movie, time_movie] = load_movie(npixels);
nimages_movie = length(time_movie);
nuclei_contours_movie = cell(nimages_movie, 1);
for i=1:nimages_movie
    B = find_active_contours(images_movie(:,:,1,i));
    nuclei_contours_movie{i} = B;
end

%%

dpERK = cat(3, squeeze(images1(:,:,2,:)),  squeeze(images2(:,:,2,:)), squeeze(images3(:,:,1,:)));
dpERK_times = [time1; time2; time3];
dpERK_contours = [nuclei_contours1; nuclei_contours2; nuclei_contours3];
% dpERK = squeeze(images3(:,:,2,:));
% dpERK_times = time3;

Dl = cat(3, squeeze(images1(:,:,3,:)),  squeeze(images2(:,:,3,:)));
Dl_times = [time1; time2];
Dl_contours = [nuclei_contours1; nuclei_contours2];

ind = squeeze(images3(:,:,2,:));
ind_times = time3;
ind_contours = nuclei_contours3;

nchannels = 3;
image_channel_sets = { dpERK; ind; Dl};
image_time_sets = {dpERK_times; ind_times; Dl_times};
image_contour_sets = { dpERK_contours; ind_contours; Dl_contours};

nuclei_channel = 1;

%%

images_movie_colored = zeros(npixels, npixels, nchannels, nimages_movie, 'uint8');
kernel_scale = 1;
weight_thres = 1e-4;
for i=1:nimages_movie
    for k=1:nchannels
        weights = exp(-(image_time_sets{k} - time_movie(i)).^2/kernel_scale.^2);
        weights = weights / sum(weights);
        
        length(weights > weight_thres)
        for j=1:length(weights)
            if weights(j) > weight_thres
                [cp1, cp2] = match_contours(image_contour_sets{k}{j}, nuclei_contours_movie{i});
                
                im1_warped = warp_image_rbf(image_channel_sets{k}(:,:,j), cp1, cp2, 0.05);
                
                images_movie_colored(:,:,k,i) = imlincomb(1, images_movie_colored(:,:,k,i), weights(j), im1_warped);
            end
        end
    end
    images_movie_colored(:,:,:,i) = imlincomb(1, images_movie_colored(:,:,:,i), 0.5, repmat(images_movie(:,:,nuclei_channel, i), 1, 1, nchannels));
end

%

% writerObj = VideoWriter('colored_movie2.avi');
% writerObj.FrameRate = 10;
% open(writerObj);
% 
% figure;
% % set(gcf, 'paperunits', 'centimeters')
% % set(gcf, 'papersize', [8 8])
% % set(gcf, 'paperposition',[0 0 8 8])
% 
% imshow(images_movie_colored(:,:,:,1), 'initialmagnification','fit','border','tight')
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
% clf
% 
% for i=1:nimages_movie
%     imshow(images_movie_colored(:,:,:,i))
%     
%     frame = getframe;
%     writeVideo(writerObj,frame);
%     clf
% end
% 
% close(writerObj);


