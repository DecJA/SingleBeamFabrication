function [maxReduction] = translation_decrease_multiSpot(xum,scaleFactor)

%Function calculates the decrease in power with distance when using the SLM
%to project multiple spots (2) using 'super'. Input is distance in meters

if nargin<2
  scaleFactor =1;
end
%Calculated on 26/09
% a1 = 0.9718;
% b1 =1.644e-05;
% c1 = 0.0002289;
%Calculated 30/09 - after slm-dist-cal
a1=12.18;
b1=-0.0005179;
c1=0.0003236;

gaussFunc = @(a,b,c,d,x) d*a1*exp(-((x-b)/c)^2);

maxReduction = gaussFunc(a1,b1,c1,scaleFactor,xum);

end