function f_rot = rotate_fn2(f, x, y, z, R)

old_size = size(f);

x = reshape(x, 1, []);
y = reshape(y, 1, []);
z = reshape(z, 1, []);
f = reshape(f, 1, []);

new_grid = R * [x; y; z];

x_rot = new_grid(1,:);
y_rot = new_grid(2,:);
z_rot = new_grid(3,:);

[theta, phi] = convert2sph(x, y, z);
[theta_rot, phi_rot] = convert2sph(x_rot, y_rot, z_rot);

f_rot = griddata(theta_rot, phi_rot, f, theta, phi, 'nearest');

f_rot = reshape(f_rot, old_size);

