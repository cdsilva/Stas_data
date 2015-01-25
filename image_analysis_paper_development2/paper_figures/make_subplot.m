function h = make_subplot(subplot_dim1, subplot_dim2, subplot_space, idx)

%[Y, X] = meshgrid((subplot_dim1:-1:1)/subplot_dim1, (1:subplot_dim2)/subplot_dim2);

X_start = mod(idx-1, subplot_dim1) / subplot_dim1 + subplot_space/2;
Y_start = 1 - ceil(idx/subplot_dim1)/subplot_dim2 + subplot_space/2;

%h = subplot('position', [X(idx)-1/subplot_dim1 Y(idx)-1/subplot_dim2 1/subplot_dim1-0.01 1/subplot_dim2-0.01]);
h = subplot('position', [X_start Y_start 1/subplot_dim1-subplot_space 1/subplot_dim2-subplot_space]);