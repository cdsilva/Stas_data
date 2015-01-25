function [corr_value, theta2_adjust] = angle_corr(theta1, theta2)

theta2_adjust = theta2 + mean(atan2d(sind(theta1-theta2), cosd(theta1-theta2)));
theta2_adjust = mod(theta2_adjust, 360);

idx = find(theta1 - theta2_adjust > 180);
theta2_adjust(idx) = theta2_adjust(idx) + 360;

idx = find(theta1 - theta2_adjust < -180);
theta2_adjust(idx) = theta2_adjust(idx) - 360;

corr_value = corr(theta1, theta2_adjust);