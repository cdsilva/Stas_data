function f_rot = rotate_fn(f, x, y, z, angle_vec)

warning('off','MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');


[x_rot, y_rot, z_rot] = rotate_grid(x, y, z, angle_vec);
%figure;
%surf(x_rot, y_rot, z_rot, f, 'linestyle', 'none')

[theta, phi] = convert2sph(x, y, z);
[theta_rot, phi_rot] = convert2sph(x_rot, y_rot, z_rot);

f_rot = griddata(theta_rot, phi_rot, f, theta, phi, 'nearest');
%figure;
%surf(x, y, z, f_rot, 'linestyle', 'cubic')
