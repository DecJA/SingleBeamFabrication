%Setup slm
addpath('R:\declan\TableD');
addpath('R:\declan\TableD\Printing');
cd('R:\declan\TableD');
run('setup_slm.m');
cd('R:\declan\TableD\Printing\CoherentProgram');
slm_check = 1;
load('zeroCorrection.mat');

if exist('d')~= 1
  %Setup DAQ board - just x/y
  d = daq.createSession('ni');
  %Make sure no voltage is output
  [chi,idi] = addAnalogInputChannel(d,'Dev2',[0:2],'Voltage');
  initialVoltage = inputSingleScan(d);
  [ch,id] = addAnalogOutputChannel(d,'Dev2',[0:2],'Voltage');
  outputSingleScan(d,[0,0,5.5]);
  %removeChannel(d,3); %Remove z-axis once calibrated
end

%Setup dithering pattern
xres = linspace(-1,1, sz(1));
freq=128;
xrescos=cos(xres.*(freq));
yrescos=cos((xres.').*freq);
xy=1*(pi/2).*sign(xrescos.*yrescos); %1 for maximum diffraction from 0 order
pattern_scatter=xy;
slm.show(pattern_scatter); 

%Get power modulation through dithering pattern
[scaleFactor,ditheredPattern] = modulateAmplitude(0.7,128,sz); 

%Draw continuous lines
ypos = 10.0e-06;
dy = 2e-6;
xlen = linspace(0,20,200).*1.0e-06;
slm.show(ditheredPattern);

yscan=1;
while yscan == 1
  y = ypos + dy;
  for ii=1:numel(xlen)
    x=xlen(ii);
    xv = dist_to_volt(x);
    yv = dist_to_volt(y);
    outputSingleScan(d,[xv,yv,initialVoltage(3)]);
    pause(0.1);
  end
end

clear all; delete(d);