clear all
close all

res = '-r300';
fmt = '-djpeg';

im_save_dir = 'paper_figures';

n = 100;

[X, Y, Z] = sphere(n);

C = zeros(size(X));
C(round(0.4*n):round(0.6*n), round(0.2*n):round(0.3*n)) = 1;

v0 = [0 0];
v1 = [-1.5 1.5];

dir = [0 0 0;
    0 1 0;
    0 0 1;
    1 0 0];

for i=1:4
    figure;
    set(gcf, 'papersize', [1 1])
    h = surf(X, Y, Z, C, 'edgecolor','none');
    set(gca,'position',[0 0 1 1],'units','normalized')
    camlight(45, -20);
    lighting phong
    material dull
    
    axis off
    hold on
    plot3(v1, v0, v0, 'k')
    plot3(v0, v1, v0, 'k')
    plot3(v0, v0, v1, 'k')
    
    if i > 1
        rotate(h,dir(i, :),25)
    end
    
    print(sprintf('%s/sphere_%d', im_save_dir, i), fmt, res)
end

figure;
for i=1:4
    figure;
    set(gcf, 'papersize', [1 1])
    h = surf(X, Y, Z, C, 'edgecolor','none');
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    camlight(45, -20);
    lighting phong
    material dull
    
    axis off
    hold on
    plot3(v1, v0, v0, 'k')
    plot3(v0, v1, v0, 'k')
    plot3(v0, v0, v1, 'k')
    if i > 1
        rotate(h,dir(i, :),25)
    else
        annotation('doublearrow', [0.4 0.6], [0.55 0.55],'color','w')
        annotation('doublearrow', [0.55 0.55], [0.4 0.6],'color','w')
        text(0, -1, 0.2, '\eta_{proj}','color','w')
        text(0.17, -1, 0, '\eta_{proj}','color','w')
        
    end
    view(0, 0)
    
    print(sprintf('%s/sphere2_%d', im_save_dir, i), fmt, res)
    
end

