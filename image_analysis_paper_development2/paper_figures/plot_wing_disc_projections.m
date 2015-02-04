function h = plot_wing_disc_projections(images, time, subplot_dim1, subplot_dim2)


nimages = size(images, 5);
npixels = size(images, 1);
nstack = size(images, 4);

for i=1:nimages
    X_start = mod(i-1, subplot_dim1) / subplot_dim1;
    Y_start = 1 - ceil(i/subplot_dim1)/subplot_dim2;

    subplot('position', [X_start Y_start 0.75/subplot_dim1 0.75/subplot_dim2]);
    imshow(recolor_wing_disc(max(images(:,:,:,:,i),[], 4)));
    if ~isempty(time)
        str1 = sprintf('%0.1f', time(i)-0.5);
        str2 = sprintf('%0.1f', time(i)+0.5);
        if strcmp(str1(end), '0');
            str1 = str1(1:end-2);
        end
        if strcmp(str2(end), '0');
            str2 = str2(1:end-2);
        end
        text(npixels, npixels, sprintf('%s-%s hr', str1, str2),'color',0.95*ones(1,3),'HorizontalAlignment','right','VerticalAlignment','bottom', 'fontsize', 6)
    end
    
    subplot('position', [X_start Y_start+0.8/subplot_dim2 0.75/subplot_dim1 0.1/subplot_dim2]);
    imshow(imresize(recolor_wing_disc(max(permute(images(:,:,:,:,i), [4 2 3 1]),[], 4)), [floor(npixels/7.5) npixels]))
    
    subplot('position', [X_start+0.8/subplot_dim1 Y_start 0.1/subplot_dim1 0.75/subplot_dim2]);
    imshow(imresize(recolor_wing_disc(max(permute(images(:,:,:,:,i), [1 4 3 2]),[], 4)), [npixels floor(npixels/7.5)]))
 
end
% h = subplot('position', [X_start Y_start 1/subplot_dim1-subplot_space 1/subplot_dim2-subplot_space]);


% 
% make_fig(17, 17/nimages);
% 
% subplot_space = 0.1;
% subplot_space2 = 0.02;
% 
% y_start = 0;
% y_mid = npixels/(npixels+nstack);
% y_end = 1;
% 
% dx1 = y_mid-y_start-subplot_space2;
% dx2 = y_end-y_mid-subplot_space;
% 
% for i=1:nimages
%     x_start = (i-1)/nimages;
%     x_mid = (i-1)/nimages+npixels/(nimages*(npixels+nstack));
%     x_end = i/nimages;
%     
%     subplot('position', [x_start y_start dx1/nimages dx1]);
%     image(recolor_wing_disc(max(images(:,:,:,:,i),[], 4)));
%     axis off
%     if ~isempty(time)
%         text(npixels/2, npixels/2, sprintf('~%2.0f hr', time(i)),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
%     end
%     
%     subplot('position', [x_start y_mid dx1/nimages dx2]);
%     image(recolor_wing_disc(max(permute(images(:,:,:,:,i), [4 2 3 1]),[], 4)))
%     axis off
%     
%     subplot('position', [x_mid y_start dx2/nimages dx1]);
%     image(recolor_wing_disc(max(permute(images(:,:,:,:,i), [1 4 3 2]),[], 4)))
%     axis off
%     
% end

function image2 = recolor_wing_disc(image)

image2 = image;

scale = 0.65;

image2(:,:,3) = immultiply(image2(:,:,3), scale);
for i=1:2
    image2(:,:,i) = imlincomb(scale, image2(:,:,i), 1, image2(:,:,3));
end
