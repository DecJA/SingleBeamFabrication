%% Main program
%Uses SLM for aberration correction and splitting the incomming beam into
%two spots. Then the stage scans the directed bitmap using the DAQ board.

slm_check = exist('sz','var');
if slm_check==0  
    addpath('R:\declan\TableD');
    addpath('R:\declan\TableD\Printing');
    cd('R:\declan\TableD');
    run('setup_slm.m');
    cd('R:\declan\TableD\Printing\CoherentProgram');
    slm_check = 1;
    load('zeroCorrection.mat');
end

%% Check if DAQ board is connected
board_check = exist('board_check','var');
calibrate_pos = 1; %Get initial voltages 
if board_check==0
  d = daq.createSession('ni');
  %Make sure no voltage is output
  [chi,idi] = addAnalogInputChannel(d,'Dev2',[0:2],'Voltage');
  initialVoltage = inputSingleScan(d);
  %removeChannels(d,1:3);
  [ch,id] = addAnalogOutputChannel(d,'Dev2',[0:2],'Voltage'); %Add 
   board_check = 1;    
end

%% Spectify grid area and print resoltion
cur=pwd;
maps=strcat(cur,'\PrintMaps');
addpath(maps);

offset =-(3/2)*pi; %Offset from SLM-Camera angle
xlength = 20e-06;
ylength = 20e-06;

%Axial steps and voltages
res_struct = 200e-09;
layers =10;

structure = 'image';%Botton left grid cornder locdddaation
im_name = 'barrier_array2';
% coords= [0e-06,2e-6];
coords = [ceil(max(xlength))/2,ceil(max(ylength)/2)];

if coords(1)<0
    disp('Pattern printing through 0 order. Position X/Y stage to be centered at 0 volts.');
elseif coords(1)>0
    disp('Pattern prining on one side of 0 order. Position X/Y stage to be at minimum at 0 volts');
end
dim = 60;
voxel_spacing = (max(xlength)/dim)*1e9;
c=newline;
disp([c,'Voxel spacing is set to: ',num2str(voxel_spacing),' nm.']);

%% Setup axial step and slm
height_struct = res_struct*layers;  
step_voltage = dist_to_volt(res_struct);

%Get current stage voltage and reposition
currentVoltage=inputSingleScan(d);
if round(currentVoltage,1)~=round(initialVoltage)
  outputSingleScan(d,initialVoltage);
end

%%  display scattering pattern on SLM to avoid polymerisation
xres = linspace(-1,1, sz(1));
freq=128;
xrescos=cos(xres.*(freq));
yrescos=cos((xres.').*freq);
xy=1*(pi/2).*sign(xrescos.*yrescos); %1 for maximum diffraction from 0 order
pattern_scatter=xy;
slm.show(pattern_scatter); 

%Sample power adjustment
reducePowerBy = 0; %0-no reduction - 1 - 100% power reduceion (works up to 0.90)
if reducePowerBy ~=0
    [scaleFactor,ditheredPattern] = modulateAmplitude(reducePowerBy,128,sz); %50% power modulation
    slm.show(ditheredPattern);
else
    ditheredPattern = zeros(sz);
end

%% Get 3d array from pattern
%Circle (cylinder), triangle, rect
%NOTE: include 'image_name','xxx' ...  if using bitmap generator
[shape_bitmap] = shape_map_2d_parse(dim,'coords',coords,...
       'xdim',xlength,'ydim',ylength,'shape','triangle');
[shape_bitmap_3d] = create_3d_structure(shape_bitmap,height_struct,res_struct);

if ceil(height_struct/res_struct) >= 2
    sz_bitmap=size(shape_bitmap_3d);
    zsteps=sz_bitmap(3);
elseif ceil(height_struct/res_struct) == 1
    zsteps=1;
end

%Display 3d bitmap to be printed
[~] = generate_3d_figure(shape_bitmap_3d);
%% Start priting pattern
xgrid=linspace(coords(1),coords(1)+xlength,dim);
ygrid=linspace(coords(2),coords(2)+ylength,dim);
zgrid=linspace(0, res_struct*layers,layers);

xplane= squeeze(shape_bitmap_3d(:,1,1));
jj=1; kk=1;

dwelTime = 0.2;
Ntotal = layers*dim*dim;

exposedPoints = length(find(squeeze(shape_bitmap_3d(:,:,1))==0));
timeLayer = (dwelTime*(dim*dim)/60)*(exposedPoints/(dim*dim));
timeEst = ((dwelTime*Ntotal)/60)*(exposedPoints/(dim*dim));

disp(['Time estimate for 1 layer: ',num2str(timeLayer), ' Min.']);
disp(['Time estimate for structure: ',num2str(timeEst), ' Min.']);

rSpot = [-10,10].*1e-6;  %static spot speration
scale_value1 = dist_to_scale(rSpot(1));
scale_value2 = dist_to_scale(rSpot(2));
sp1 = otslm.simple.linear(sz,scale_value1,'angle',theta1+offset);
sp2 = otslm.simple.linear(sz,scale_value2,'angle',theta2+offset);
pattern = otslm.tools.combine({sp1,sp2},'method','super');

%To Do: include distance-power modulation
pattern = otslm.tools.finalize(pattern);
allPatterns = pattern + zeroCorrection + ditheredPattern;
slm.show(allPatterns);

if coords(1) > 0
    structLength = (abs(min(xlength) - max(xlength)));
elseif coords(1)<0
    structLength = (abs(min(xlength) - max(xlength)))/2;
end

if structLength>(abs(rSpot(1)-rSpot(2)))
    warning('Cancel Printing: Structures will overlap')
else
    structSep = structLength - abs(rSpot(1)-rSpot(2));
    disp(['Structure seperation: ',num2str(structSep.*1e6),' Microns.']);
    return
end
    
for ii=1:layers
  z = zgrid(ii);
  for jj=1:dim
    y = ygrid(jj);
    for kk=1:dim
      
      x=xgrid(kk);
      [xv,yv,zv] = dist_to_volt_3d(x,y,z);
      pos = [xv,yv,-zv];      
      
      if shape_bitmap_3d(kk,jj,ii)==1
        outputSingleScan(d,pos);
        pause(dwelTime);
      end 
      %add: zero order correcting, voltage tracking, motion sensor,
      %abberration correction
    end
  end
     disp(['Moving to step: ',num2str(ii),' of', num2str(zsteps)]);
end
slm.show(pattern_scatter);
disp(['Program Completed.']);

