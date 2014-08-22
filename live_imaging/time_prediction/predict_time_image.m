function [t, tmin, tmax] = predict_time_image(image)
% this function takes an image and returns the predicted time, along with
% the 95% confidence intervals for the time
% image is a single-channel image of nuclei, oriented with the dorsal side
% at the top
% t is the predicted time
% tmin is the minimum predicted time in the 95% confidence interval
% tmax is the maximum predicted time in the 95% confidence interval

persistent mean_train_data V_train mdl;  % Declare data as a persistent variable
if isempty(mean_train_data) || isempty(V_train) || isempty(mdl)
    load('PCA_train_data.mat');
end

rot_angle = 0;
data = (create_PCA_data(image_normalize(image, rot_angle)) - mean_train_data)*V_train;

[ypred,yci] = predict(mdl,data, 'prediction','observation');

t = ypred;
tmin = yci(1);
tmax = yci(2);


