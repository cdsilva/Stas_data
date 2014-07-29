clear all
close all

clear all
close all

%%

% image_fn = @(image) medfilt2(adapthisteq(crop_image(image, 200, 0)),[5 5]);
image_fn = @(image) adapthisteq(crop_image(image, 200, 0));
% image_fn = @(image) adapthisteq(image);

channel = 1;
npixels = 100;

movie1 = '14_0624/emb01_hisRFP_gastrulation.avi';
movie2 = 'bomyi_emb01_gast01.avi';

%%
[images1, time1] = read_video(movie1, npixels);
images1 = images1(:, :, channel, :);
images1 = squeeze(images1);

nimages1 = length(time1);
subplot_dim1 = ceil(sqrt(nimages1));
subplot_dim2 = ceil(nimages1 / subplot_dim1);


[images2, time2] = read_video(movie2, npixels);
images2 = images2(:, :, channel, :);
images2 = squeeze(images2);

nimages2 = length(time2);
subplot_dim1_movie = ceil(sqrt(nimages2));
subplot_dim2_movie = ceil(nimages2 / subplot_dim1_movie);

%%

for i=1:nimages1
    images1(:,:,i) = image_fn(images1(:,:,i));
end

for i=1:nimages2
    images2(:,:,i) = image_fn(images2(:,:,i));
end

figure;
for i=1:nimages1
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(images1(:,:,i));
end


figure;
for i=1:nimages1
    subplot(subplot_dim1, subplot_dim2, i)
    imshow(adapthisteq(images1(:,:,i)-uint8(0.5*mean(double(images2), 3))));
end

figure;
imshow(uint8(mean(double(images2), 3)))

figure;
im_tmp = zeros(npixels, npixels, 3, 'uint8');
im_tmp(:,:,1) = images1(:,:,30);
im_tmp(:,:,2) = images2(:,:,30);
imshow(im_tmp)

%%

data1 = reshape(double(images1), npixels^2, [])';
data2 = reshape(double(images2), npixels^2, [])';

mean1 = mean(data1);
mean2 = mean(data2);

data1 = data1 - repmat(mean1, nimages1, 1);
data2 = data2 - repmat(mean2, nimages2, 1);

%%

% [V, D, coeff] = PCA_images(data1, 100);
[u, s, V, eigenimages1, D, coeff] = PCA_images(images1, nimages1);

figure;
bar(diag(D) / sum(diag(D)));

figure;
plot(cumsum(diag(D)) / sum(diag(D)),'.');

figure;
for i=1:4
    subplot(2,2,i)
    imshow(eigenimages1(:,:,i))
end

nmodes = 23;
images1_denoised = zeros(size(images1), 'uint8');
for i=1:nimages1
    %     images1_denoised(:,:,i) = reshape(uint8(u(i, 1:nmodes)*s(1:nmodes,1:nmodes)*V(:, 1:nmodes)'+mean1), [npixels npixels]);
    images1_denoised(:,:,i) = reshape(uint8(data1(i,:)*V(:, 1:nmodes)*V(:, 1:nmodes)'+mean1), [npixels npixels]);
end

images2_denoised = zeros(size(images2), 'uint8');
for i=1:nimages2
    images2_denoised(:,:,i) = reshape(uint8(data2(i,:)*V(:, 1:nmodes)*V(:, 1:nmodes)'+mean2), [npixels npixels]);
end

%%

figure;
writerObj = VideoWriter('movie1_pca.avi');
writerObj.FrameRate = 15;
open(writerObj);
imshow(images1_denoised(:,:,1));
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
for i=1:nimages1
    imshow(images1_denoised(:,:,i));
   frame = getframe;
   writeVideo(writerObj,frame);
end
close(writerObj);

figure;
writerObj = VideoWriter('movie2_pca.avi');
writerObj.FrameRate = 15;
open(writerObj);
imshow(images2_denoised(:,:,1));
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
for i=1:nimages2
    imshow(images2_denoised(:,:,i));
   frame = getframe;
   writeVideo(writerObj,frame);
end
close(writerObj);

%%

figure;
for i=1:4
    subplot(2,2,i)
    plot(data1*V(:, i),'.')
    xlabel('time')
    ylabel(sprintf('data 1 onto mode %d', i))
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(data2*V(:, i),'.')
    xlabel('time')
    ylabel(sprintf('data 2 onto mode %d', i))
end


