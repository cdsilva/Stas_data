function sx = compute_scattering_coeff(image)

[npixels, ~] = size(image);

% set scattering transform parameters
filt_opt = struct();
filt_opt.J = ceil(log2(npixels));
filt_opt.L = 8;

filt_rot_opt = struct();

scat_opt = struct();
% scat_opt.oversampling = 100;

% compute scattering transform of first image
x = double(image);
% define wavelet transforms
Wop = wavelet_factory_3d(size(x), filt_opt, filt_rot_opt, scat_opt);

Sx = scat(x, Wop);
sx = mean(mean(format_scat(Sx),2),3)';
% sx = mean(format_scat(Sx),3);
% sx = format_scat(Sx);
% sx = sx(:);
