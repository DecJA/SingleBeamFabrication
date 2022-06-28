function [xv,yv,zv] = dist_to_volt_3d(x,y,z)

%Converts requested chnage in distance to DAQ acquisition card voltage
%change. Moves in relative distance only, not absolute. Voltage offset does
%not change with manual DC offset on controller

%Madcity labs 3 axis stage
mx = 3e-5; %Experimentally verified slopes
my = mx;
mz = mx;

xv = x/mx;
yv = y/my;
zv = z/mz;

end