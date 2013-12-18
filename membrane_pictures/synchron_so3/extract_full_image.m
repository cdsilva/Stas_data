function im1 = extract_full_image(f, im_size)

scale = 4;

nim1 = floor(im_size(1)/(scale-2))*(scale-2);
nim2 = floor(im_size(2)/(scale-2))*(scale-2);

npoints1 = nim1/(scale-2);
npoints2a = nim2/(scale-2);
npoints2b = nim2*(scale+1)/(scale-2);

im1 = f(npoints1+1:npoints1+nim1,npoints2a+1:npoints2a+nim2);