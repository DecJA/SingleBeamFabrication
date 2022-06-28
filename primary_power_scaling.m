function [scaleFactor] = primary_power_scaling(P0,x)

% a1 = 0.9718;
% b1 =1.644e-05;
% c1 = 0.0002289;
%Calculated 30/09 - after slm-dist-cal
a1=12.18;
b1=-0.0005179;
c1=0.0003236;

scaleFactor = (P0/a1)*exp(((x-b1)/c1)^2);

end