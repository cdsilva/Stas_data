function [theta, phi] = convert2sph(x, y, z)

theta = acos(z);
phi = atan2(y, x);