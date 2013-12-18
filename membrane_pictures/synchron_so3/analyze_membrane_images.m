% read in membrane pictures; align using vector dmaps

clear all
close all

res = '-r300';
fmt = '-djpeg';

cm_green = gray;
cm_green(:,1) = 0;
cm_green(:,3) = 0;

% set parallel parameters
run_in_parallel = false;
nprocs = 4;

%set image parameters
nimages = 25;
subplot_dim1 = 5;
subplot_dim2 = 5;

image_dir = '../membrane2';
adjust_image = true;
image_channel = 2;
im_scale = 0.2;
max_deg = 7;
% load in membrane length data
load('../membrane2/length_oct16.mat');
mem_lengths = length(:,2);
clear length;

% dimension of rotations
rot_dim = 3;

% where to store fourier expansion coefficients
ncoeff = sum(2*(0:max_deg)+1);
F_l = zeros(nimages, ncoeff);

image_set = cell(nimages,1);

% load in images
figure;
for i=1:nimages
    im1 = imread(sprintf('%s/emb%02d.tif', image_dir, i));
    im1 = imresize(im1, im_scale);
    im1 = im1(:,:,image_channel);
    if adjust_image
        im1 = imadjust(im1);
        im1 = adjust_image_radially(im1);
        %im1 = medfilt2(im1, [3 3]);
        %im1 = edge(im1);
        %h = fspecial('average',[2 2]);
        %im1 = imfilter(im1,h);
        %im1 = im2bw(im1, 0.25);
    end
    subplot(subplot_dim1,subplot_dim2,i)
    imshow(im1)
    im1 = double(im1);
    image_set{i} = im1;
    
    temp_F_l = zeros(1, ncoeff);
    
    for l=0:max_deg
        start_idx = sum(2*(0:(l-1))+1)+1;
        temp_F_l(start_idx:start_idx+2*l) = f_hat(l, im1);
    end
    
    F_l(i,:) = temp_F_l;
end

%% plot rotation of images
n = size(im1, 1);
[x, y, z] = sphere(2*n);
theta = acos(z);
phi = atan2(y, x);
f = zeros(size(theta));
%
% test_ind = 4;
% for l=0:max_deg
%     start_idx = sum(2*(0:(l-1))+1)+1;
%     sph_harm = Y_lm(l, reshape(theta, 1, []), reshape(phi, 1, []));
%     f = f + reshape(sum(repmat(F_l(test_ind,start_idx:start_idx+2*l).', 1, size(sph_harm,2)).*conj(sph_harm), 1), size(f));
% end
%
% figure;
% surf(x, y, z, real(f), 'linestyle','none')
% colormap(gray)
%
% [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, [0 pi/4 0]);
% figure;
% surf(x_rot, y_rot, z_rot, real(f), 'linestyle','none')
% colormap(gray)
%
% [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, [0 pi/2 0]);
% figure;
% surf(x_rot, y_rot, z_rot, real(f), 'linestyle','none')
% colormap(gray)
%
% [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, [0 3*pi/4 0]);
% figure;
% surf(x_rot, y_rot, z_rot, real(f), 'linestyle','none')
% colormap(gray)
%
% return

%% compute pairwise alignments

H = zeros(rot_dim*nimages);
H2 = zeros(nimages, rot_dim*nimages, rot_dim);

W = zeros(nimages);
opt_angles = zeros(nimages, nimages, rot_dim);

if run_in_parallel
    matlabpool('open',nprocs);
end

parfor i=1:nimages
    i
    temp_opt_angles = zeros(1, nimages, rot_dim);
    temp_H = zeros(rot_dim, rot_dim*nimages);
    temp_W = zeros(1, nimages);
    for j=1:i-1
        [F_l_aligned, A] = align_fhat(F_l(j,:), F_l(i,:), max_deg);
        
        temp_opt_angles(1, j, :) = A;
        
        R = rot_matrix(A);
        
        temp_H(:, rot_dim*j-(rot_dim-1):rot_dim*j) = R;
        
        %temp_W(j) = norm(F_l_aligned-F_l(i,:))^2;
        im1 = image_set{i}
        im2 = image_set{j};

        [x1, y1, z1, f1] = project_full_image(im1);
        [x2, y2, z2, f2] = project_full_image(im2);
        f2 = rotate_fn(f2, x2, y2, z2, A);
        temp_W(j) = norm(f1-f2).^2
    end
    
    opt_angles(i,:,:) = temp_opt_angles;
    W(i,:) = temp_W;
    for k=1:rot_dim
        H2(i,:,k) = temp_H(k,:);
    end
    
end

for k=1:rot_dim
    H(k:rot_dim:end,:) = H2(:,:,k);
end
opt_angles = opt_angles - permute(opt_angles, [2 1 3]);
H = H + H';
W = W + W';

if run_in_parallel
    matlabpool close;
end

%save('opt_alignments.mat','H','W','opt_angles','F_l','mem_lengths','nimages','im_scale','max_deg');
%return;


%% angular synchronization
[V, D] = eigs(H, 3);

% find optimal rotations
[u, s, v] = svd(V(1:3,1:3));
R0 = u*v';
R = R0;
%R0 = [-1 0 0; 0 -1 0; 0 0 1]*R0;
figure;
for i=1:nimages
    [u, s, v] = svd(V(3*i-2:3*i,1:3));
    R = u*v';
    if i==1
        angle_vec = get_angles(eye(3));
    else
        angle_vec = get_angles(R0*R');
    end
    
    im1 = image_set{i};
    [x2, y2, z2, f2] = project_full_image(im1);
    f_rot = rotate_fn(f2, x2, y2, z2, angle_vec);
    %[x_rot, y_rot, z_rot] = rotate_grid(x2, y2, z2, R0*R');
    
    subplot(subplot_dim1,subplot_dim2,i)
    surf(x2, y2, z2, f_rot, 'linestyle','none')
    colormap(gray)
    view(180,0)
    axis equal
end

Waligned = zeros(nimages);
for i=1:nimages
    [u, s, v] = svd(V(3*i-2:3*i,1:3));
    R = u*v';
    if i==1
        angle_vec = get_angles(eye(3));
    else
        angle_vec = get_angles(R0*R');
    end
    
    im1 = image_set{i};
    [x2, y2, z2, f] = project_full_image(im1);
    f_rot1 = rotate_fn(f, x2, y2, z2, angle_vec);
    for j=1:i
        [u, s, v] = svd(V(3*j-2:3*j,1:3));
        R = u*v';
        if j==1
            angle_vec = get_angles(eye(3));
        else
            angle_vec = get_angles(R0*R');
        end
        
        im1 = image_set{j};
        [x2, y2, z2, f] = project_full_image(im1);
        f_rot2 = rotate_fn(f2, x2, y2, z2, angle_vec);
        Waligned(i,j) = norm(f_rot1-f_rot2)^2;
        Waligned(i,j) = Waligned(j,i);
    end
end


%% vector DMAPS
eps = 5e7;
W1 = exp(-W/eps);
W1 = diag(1./sum(W1)) * W1;
A = zeros(size(H));
for i=1:nimages
    for j=1:nimages
        A(3*i-2:3*i,3*j-2:3*j) = H(3*i-2:3*i,3*j-2:3*j) * W1(i,j);
    end
end

neigs = 5;
[V2, D2] = eigs(A, neigs);

%% find optimal rotations
% [u, s, v] = svd(V(1:3,1:3));
% R0 = u*v';
% R0 = [-1 0 0; 0 -1 0; 0 0 1]*R0;
% for i=1:nimages
%     f1 = zeros(size(theta));
%     for l=0:max_deg
%         start_idx = sum(2*(0:(l-1))+1)+1;
%         sph_harm = Y_lm(l, reshape(theta, 1, []), reshape(phi, 1, []));
%         f1 = f1 + reshape(sum(repmat(F_l(i,start_idx:start_idx+2*l).', 1, size(sph_harm,2)).*conj(sph_harm), 1), size(f1));
%     end
%
%     [u, s, v] = svd(V(3*i-2:3*i,1:3));
%     R = u*v';
%     [x_rot, y_rot, z_rot] = rotate_grid(x, y, z, R0*R');
%     figure;
%     surf(x_rot, y_rot, z_rot, real(f1), 'linestyle','none')
%     colormap(gray)
%
%     im1 = imread(sprintf('%s/emb%02d.tif', image_dir, i));
%     im1 = imresize(im1, im_scale);
%     im1 = im1(:,:,image_channel);
%     if adjust_image
%         im1 = imadjust(im1);
%     end
%     im1 = double(im1);
%     [x2, y2, z2, f2] = project_full_image(im1);
%     [x_rot, y_rot, z_rot] = rotate_grid(x2, y2, z2, R0*R');
%     figure;
%     surf(x_rot, y_rot, z_rot, f2, 'linestyle','none')
%     colormap(gray)
% end

%% find embedding

V_coord = zeros(nimages, neigs*(neigs-1)/2);

for i=1:nimages
    coord_idx = 1;
    for j1=1:neigs
        for j2=j1:neigs
            V_coord(i,coord_idx) = V2(3*i-2:3*i,j1)'*V2(3*i-2:3*i,j2);
            coord_idx = coord_idx + 1;
        end
    end
end

for i=1:size(V_coord,2)
    figure;
    plot(V_coord(:,i),mem_lengths(1:nimages),'.')
end

V_coord_idx = 4;
[~, sort_idx] = sort(V_coord(:,V_coord_idx));
%select_images = [1 3 6 9 14];
select_images = [1:5];
figure;
for j=1:5
    i = sort_idx(select_images(j));
    [u, s, v] = svd(V(3*i-2:3*i,1:3));
    R = u*v';
    
    im1 = image_set{i};
    [x2, y2, z2, f2] = project_full_image(im1);
    [x_rot, y_rot, z_rot] = rotate_grid(x2, y2, z2, R0*R');
    
    subplot(4,5,j)
    surf(x_rot, y_rot, z_rot, f2, 'linestyle','none')
    colormap(cm_green)
    view(180,0)
    axis equal
end
subplot(4,5,6:20)
plot(V_coord(:,V_coord_idx),mem_lengths(1:nimages),'.')
hold on
plot(V_coord(sort_idx(select_images),V_coord_idx),mem_lengths(sort_idx(select_images)),'or')
xlabel('embedding coordinate from VDM')
ylabel('membrane length measured from experiments')
%print('membrane_images_synch2',fmt,res)



