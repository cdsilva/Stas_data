clear all
close all

res = '-r300';
fmt = '-djpeg';

%set(0,'DefaultLineMarkerSize',8)
%set(0,'DefaultLineLineWidth',2)
%set(0,'DefaultAxesFontSize',20)

im_save_dir = 'paper_figures2';

n = 25;

[X, Y, Z] = sphere(n);

C = zeros(size(X));
C(round(0.4*n):round(0.6*n), round(0.2*n):round(0.3*n)) = 1;
%C = Z;
%C(round(0.4*n):round(0.6*n), round(0.2*n):round(0.3*n)) = -2;

v0 = [0 0];
v1 = [-1.5 1.5];

dir = [0 0 0;
    0 1 0;
    0 0 1;
    1 0 0];
theta = 0:0.1:pi;
r = 0.2;
r1 = 1.4;

for i=1:4
    figure;
    set(gcf, 'units', 'centimeters')
    set(gcf, 'papersize', [4 4])
    set(gcf, 'paperposition', [0 0 4 4])
    h = surf(X, Y, Z, C, 'edgecolor','k');
     set(gca,'position',[0 0 1 1],'units','normalized')
    camlight(45, -20);
    lighting phong
    %material dull
    set(h, 'diffusestrength', 1);
    set(h, 'specularstrength', 1);
    %set(gca, 'ambientlightcolor', [1 0 0]);
    set(h, 'AmbientStrength', 0.7);
    %set(h, 'SpecularColorReflectance', 0.25);
    %set(h, 'SpecularExponent', 5.0);
    
    axis off
    hold on
    plot3(v1, v0, v0, 'k', 'linewidth', 2)
    plot3(v0, v1, v0, 'k', 'linewidth', 2)
    plot3(v0, v0, v1, 'k', 'linewidth', 2)
    
    if i > 1
        rotate(h,dir(i, :),25)
    end
    if i == 2
        plot3(r*cos(theta), -1.75*ones(size(theta)),r*sin(theta), 'linewidth', 2)
        plot3(r*cos(theta(1)), -1.75,r*sin(theta(1)), '<', 'linewidth', 2, 'markerfacecolor','b')
        text(r+0.1, -r1,0.1,  '\alpha', 'color','b', 'fontsize', 24)
            plot3(v0, 1.1*v1, v0, 'k', 'linewidth', 2)
    axis equal
    end
    
    if i == 3
        plot3(r*cos(theta), r*sin(theta), -r1*ones(size(theta)), 'linewidth', 2)
        plot3(r*cos(theta(end)), r*sin(theta(end)), -r1, '<', 'linewidth', 2, 'markerfacecolor','b')
        text(0, r+0.02, -r1-0.2, '\beta', 'color','b', 'fontsize', 24)
    end
    if i == 4
        plot3(-r1*ones(size(theta)), r*sin(theta+pi), r*cos(theta+pi), 'linewidth', 2)
        plot3(-r1, r*sin(theta(1)+pi), r*cos(theta(1)+pi), '<', 'linewidth', 2, 'markerfacecolor','b')
        text(-r1, -r, -r-0.1, '\gamma', 'color','b', 'fontsize', 24)
    end
    
zoom(1.5)
print(sprintf('%s/sphere_%d', im_save_dir, i), fmt, res)
end

theta = pi/4:0.1:3*pi/2;
arrow_slope = 600;
r1 = 1.25;

for i=1:4
    figure;
    set(gcf, 'units', 'centimeters')
    set(gcf, 'papersize', [4 4])
    set(gcf, 'paperposition', [0 0 4 4])
    h = surf(X, Y, Z, C, 'edgecolor','k');
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    camlight(-10, 30);
    lighting phong
    %material dull
    set(h, 'diffusestrength', 3.0);
    set(h, 'specularstrength', 1);
    %set(gca, 'ambientlightcolor', [1 0 0]);
    set(h, 'AmbientStrength', 0.7);
    %set(h, 'SpecularColorReflectance', 0.25);
    %set(h, 'SpecularExponent', 5.0);
    
    axis off
    hold on
    plot3(v1, v0, v0, 'k')
    plot3(v0, v1, v0, 'k')
    plot3(v0, v0, v1, 'k')
    if i > 1
        rotate(h,dir(i, :),25)
    else
        %annotation('doublearrow', [0.4 0.6], [0.55 0.55],'color','w')
        %annotation('doublearrow', [0.55 0.55], [0.4 0.6],'color','w')
        plot3([-0.53 0.42], [-2 -2], [0.1 0.1], 'color',0.8*ones(1,3), 'linewidth', 2)
        plot3([0.1 0.1], [-2 -2], [-0.41 0.31], 'color',0.8*ones(1,3), 'linewidth', 2)
        plot3(-0.51, -2, 0.1, '<','color',0.8*ones(1,3), 'linewidth', 2, 'markerfacecolor',0.8*ones(1,3))
        plot3(0.4, -2, 0.1, '>','color',0.8*ones(1,3), 'linewidth', 2, 'markerfacecolor',0.8*ones(1,3))
        plot3(0.1, -2, -0.39, 'v','color',0.8*ones(1,3), 'linewidth', 2, 'markerfacecolor',0.8*ones(1,3))
        plot3(0.1, -2, 0.3, '^','color',0.8*ones(1,3), 'linewidth', 2, 'markerfacecolor',0.8*ones(1,3))
        text(-0.5, -1, 0.25, '\eta_{proj}','color',0.8*ones(1,3), 'fontsize', 20)
        text(0.17, -1, -0.2, '\eta_{proj}','color',0.8*ones(1,3), 'fontsize', 20)
    end
    if i == 2
        plot3(r*cos(theta), -r1*ones(size(theta)),r*sin(theta), 'linewidth', 2, 'color', 0.8*ones(1,3))
        plot3(r*cos(theta(1)), -r1,r*sin(theta(1)), '<', 'linewidth', 2, 'color', 0.8*ones(1,3),'markerfacecolor',0.8*ones(1,3))
        text(-r-0.1, -r1, r+0.1,  '\alpha', 'color',0.8*ones(1,3), 'fontsize', 28)
        
    end
    
    if i == 3
        plot3(r*cos(theta), r*sin(theta), -r1+(1:length(theta))/arrow_slope, 'linewidth', 2)
        plot3(r*cos(theta(end)), r*sin(theta(end)), -r1+length(theta)/arrow_slope, '>', 'linewidth', 2, 'markerfacecolor','b')
        text(0, r+0.05, -r1-0.1, '\beta', 'color','b', 'fontsize', 24)
    end
    if i == 4
        plot3(-r1+(1:length(theta))/arrow_slope, r*sin(theta+pi), r*cos(theta+pi), 'linewidth', 2)
        plot3(-r1+length(theta)/arrow_slope, r*sin(theta(end)+pi), r*cos(theta(end)+pi), 'v', 'linewidth', 2, 'markerfacecolor','b')
        text(-r1, -r, -r-0.1, '\gamma', 'color','b', 'fontsize', 24)
    end
    view(0, 0)
    
    print(sprintf('%s/sphere2_%d', im_save_dir, i), fmt, res)
    
end

