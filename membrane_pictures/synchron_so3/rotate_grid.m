function [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, angles)

if size(angles,1) == 1 || size(angles,2)==1
    R = rot_matrix(angles);
else
    R = angles;
end

points = [reshape(x, 1, []); reshape(y, 1, []); reshape(z, 1, [])];
points = R*points;
x_rot = reshape(points(1,:),size(x));
y_rot = reshape(points(2,:),size(y));
z_rot = reshape(points(3,:),size(z));