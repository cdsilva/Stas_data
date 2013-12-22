function y = power_law_params(time, coeff)
% this function finds the optimal y to fit coeff = y(1) * time^y(2)

options = optimset('display','off');
%y = lsqcurvefit(@power_law, [0.1; 3; -400], time, coeff, [], [], options);
y = lsqcurvefit(@power_law, [1e-4; 3], time, coeff, [], [], options);