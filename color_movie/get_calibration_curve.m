clear all
close all

A = imread('Ingression_kinetics_Anna Sokac.tif');


xlim = [0 80];
ylim = [0 60];

figure; 
imshow(A(128:669,161:850,:), 'xdata', xlim, 'ydata',ylim )

[x,y] = ginput;

cal_curve_time = x;
cal_curve_membrane = y;

save('calibration_curve_data', 'cal_curve_time','cal_curve_membrane');

