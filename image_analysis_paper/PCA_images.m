function [eigenimages, D, proj_coeff] = PCA_images(image_set, neigenimages)

[npixels, npixels2, nchannels, m] = size(image_set);

data = zeros(m, npixels*npixels2*nchannels);

for i=1:m
    im1 = image_set(:,:,:,i);
    data(i, :) = double(im1(:)');
end

data_mean = mean(data);
data = data - repmat(data_mean, m, 1);
%[V, D] = PCA(data, neigenimages);
[u, s, V] = svds(data, neigenimages);

D = s.^2;

max_intensity = 255;
eigenimages = zeros(npixels, npixels2, nchannels, neigenimages, 'uint8');

for i=1:neigenimages
    im1 = V(:,i);
    if max(im1) < -min(im1)
        im1 = -im1;
    end
    %im1 = im1 / max(im1) * (max_intensity - max(data_mean));
    %im1 = im1 + data_mean';
    im1 = abs(im1 / max(im1) * max_intensity);
    eigenimages(:, :, :, i) = uint8(reshape(im1, npixels2, npixels2, nchannels));
end

%proj_coeff = data * V;
proj_coeff = u*s;

