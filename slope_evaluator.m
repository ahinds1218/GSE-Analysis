function [y_smooth,y_smooth_prime] = slope_evaluator( x,y )
% slope of {x,y} using csaps cubic smoothing spline.
% based on TEST_csapi.m on Olin 213 lab computer Lenovo M92
% folder E:\Ud\Matlab\Adam\Drive8_code
% _2013_08_09
% _2013_08_10
% -2013_08_13

% inputs: x(j), y(:,j) as required by csapi 
w = linspace(0.2,0.2,length(x));  % weights.
w(1)   = 0.05;
w(end) = 0.05;
w(7) = .05;
w(8) = .05;
xx = linspace(0,1.1*x(end),151);

% --apply cubic spline smoothing function csaps
[pp2, psmooth ] = csaps(x,y);  % use csaps to get p
p_new = 0.9999*psmooth;        % adjust p slightly
pp3 = csaps(x,y,p_new,[],w);  % use csaps again to get pp
y_smooth = (fnval(pp3,x)); % pp3 = csaps(x,y,p_new,[],w);

%--- evaluate first derivative of smoothed data
p_der_3= fnder(pp3,1);
y_smooth_prime = (ppval(p_der_3,x));

end

