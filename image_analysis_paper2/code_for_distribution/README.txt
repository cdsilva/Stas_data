This directory contains matlab code for image registration and temporal ordering, 
as described in 
"Temporal ordering and registration of images in studies of developmental dynamics"

All code is implemented in MATLAB, and was tested using MATLAB R2013b. 

There are two example scripts:
analyze_live_imaging.m reads in a live imaging movie (also provided), and registers and orders the frames
analyze_fixed_images.m reads in a set of fixed images (also provided), registers and orders them, and also computes an average trajectory of these images. 

To register and order a set of images, the workflow is as follows:
1. Read in and store the images in an npixels x npixels x nchannels x nimages array (for multichannel images), 
or an npixels x npixels x nimages array (for grayscale images). 
All images should be preprocessed as necessary 
(signals whose absolute value is not informative should be normalized, images should be resized/cropped if size differences are not important, etc.)

2. Call the "compute_pairwise_alignments" function to compute the pairwise alignments (R) and distances (W) between images. 

3. Call the "vdm" function, using the rotations R and the distance matrix W computed in the previous step, 
to compute the optimal rotations (R_opt) and the embedding coordinates (embed_coord) using vector diffusion maps.

4. Call the "register_all_images" function, using the optimal rotations R_opt computed by vdm in the previous step, to register the images. 

5. Order the registerd images by sorting them by the value of the first embedding coordinate, embed_coord(:,1).

6. If desired, compute an average trajectory by calling the "compute_average_trajectory" function, using the set of registered and ordered images as input. 

A few notes: 
1. The images are registered up to a constant rotation (and translation, if relevant).
After registration, the entire set can be rotated by a constant angle if desired. 

2. The sign of the embedding coordinates is not meaningful and can be switched.
This corresponds to reversing the order of the images. 

3. Higher image resolution will make the code slower, 
and so the image resolution should be chosen as low as possible, while still retaining the important image features.

4. Including translations also increases the computational cost. 
