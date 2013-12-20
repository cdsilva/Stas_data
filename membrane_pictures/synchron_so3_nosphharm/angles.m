function [alpha, beta, gamma] = angles(R)

alpha = atan2(R(3,2), R(3,3));
gamma = atan2(R(2,1), R(1,1));
beta = -asin(R(3,1));
