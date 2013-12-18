
res = '-r300';
fmt = '-djpeg';

cm_red = gray;
cm_red(:,2) = 0;
cm_red(:,3) = 0;

% find optimal rotations
[u, s, v] = svd(V(1:3,1:3));
R0 = u*v';
%R0 = [-1 0 0; 0 -1 0; 0 0 1]*R0;
fig1 = figure;
fig2 = figure;
for i=1:nimages
    [u, s, v] = svd(V(3*i-2:3*i,1:3));
    R = u*v';
    [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, R0*R');
    
    im1 = image_set{i};
    [x2, y2, z2, f2] = project_full_image(im1);
    figure;
    surf(x2, y2, z2, f2, 'linestyle','none')
    colormap(cm_red)
    view(180,0)
    axis equal
    axis off
    print(sprintf('dpERK_unaligned_%d',i),fmt,res)

    [x_rot, y_rot, z_rot] = rotate_grid(x2, y2, z2, R0*R');
    figure;
    surf(x_rot, y_rot, z_rot, f2, 'linestyle','none')
    colormap(cm_red)
    view(180,0)
    axis equal
    axis off
    print(sprintf('dpERK_aligned_%d',i),fmt,res)
end