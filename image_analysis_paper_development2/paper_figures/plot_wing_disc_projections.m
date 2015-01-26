function h = plot_wing_disc_projections(images, time)


nimages = size(images, 5);
npixels = size(images, 1);
nstack = size(images, 4);

make_fig(17, 17/nimages);

subplot_space = 0.1;
subplot_space2 = 0.02;

y_start = 0;
y_mid = npixels/(npixels+nstack);
y_end = 1;

dx1 = y_mid-y_start-subplot_space2;
dx2 = y_end-y_mid-subplot_space;

for i=1:nimages
    x_start = (i-1)/nimages;
    x_mid = (i-1)/nimages+npixels/(nimages*(npixels+nstack));
    x_end = i/nimages;
    
    
    
    %     nimages*(x_mid-x_start-subplot_space2)
    %     y_mid-y_start-subplot_space2
    
    subplot('position', [x_start y_start dx1/nimages dx1]);
    image(max(images(:,:,:,:,i),[], 4));
    axis off
    if ~isempty(time)
        text(npixels/2, npixels/2, sprintf('~%2.0f hr', time(i)),'color',0.95*ones(1,3),'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 6)
    end
    
    subplot('position', [x_start y_mid dx1/nimages dx2]);
    image(max(permute(images(:,:,:,:,i), [4 2 3 1]),[], 4))
    axis off
    
    subplot('position', [x_mid y_start dx2/nimages dx1]);
    image(max(permute(images(:,:,:,:,i), [1 4 3 2]),[], 4))
    axis off
    
end
