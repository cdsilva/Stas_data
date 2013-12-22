function x = inv_power_law(params, y)
% this function returns the value of the inverse power law at the point y

x = (y / params(1)).^(1/params(2));