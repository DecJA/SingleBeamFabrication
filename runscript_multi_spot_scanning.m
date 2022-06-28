
slm_check = exist('sz','var');
if slm_check==0  
    addpath('R:\declan\TableD');
    addpath('R:\declan\TableD\Printing');
    cd('R:\declan\TableD');
    run('setup_slm.m');
    cd('R:\declan\TableD\Printing\CoherentProgram');
    slm_check = 1;
    offset = -(3/2)*pi;
    load('zeroCorrection_191121.mat');
%     load('zeroCorrection.mat');
end

%% Check if DAQ board is connected
board_check = exist('board_check','var');
startVoltage = [0,0,0];
if board_check==0
  d = daq.createSession('ni');
  %Make sure no voltage is output
  [chi,idi] = addAnalogInputChannel(d,'Dev2',[0:2],'Voltage'); %3-axis
  [ch,id] = addAnalogOutputChannel(d,'Dev2',[0:2],'Voltage'); %3-axis
  outputSingleScan(d,startVoltage); %Start with large positive voltage for axial dim.
  initialVoltage = inputSingleScan(d);
  board_check = 1;  
end

%% Spectify grid area and print resoltion
cur=pwd;
maps=strcat(cur,'\PrintMaps');
addpath(maps);

offset =-(3/2)*pi; %Offset from SLM-Camera angle
xlength = 5e-06;
ylength = 5e-06;

%Axial steps and voltages
dim =40;
layers = 10;
res_struct = 200e-09;
dPz=(linspace(0,.2.^1,layers+1));
dPz(1)=-.2;
dPz(2)=-.1;
dPz(7:end)=dPz(7);
dwelTime =0.04;
reducePowerBy = 0.1; %0-no reduction - 1 - 100% power reduceion (works up to 0.90)


structure = 'image';%Botton left grid cornder locdddaation
im_name = 'flower_bulk2';
coords= [2e-06,0e-06];
voxel_spacing = (max(xlength)/dim)*1e9;
c=newline;
disp([c,'Voxel spacing is set to: ',num2str(voxel_spacing),' nm.']);

%Multiple spot positions
rSpot = [7,7,7,7].*1e-6;
thetaSpot = [2*pi,pi/2,-1*pi,-pi/2];
%% Setup axial step and slm
height_struct = res_struct*layers;  
step_voltage = dist_to_volt(res_struct);

maxAxialVoltage = (res_struct*step_voltage);
if maxAxialVoltage >= (10)
  warning('WARNING: Total axial displacement exceedes travel limit.');
end
%%  display scattering pattern on SLM to avoid polymerisation
%   calculate dithering pattern to scale power
xres = linspace(-1,1, sz(1));
freq=256;
xrescos=cos(xres.*(freq));
yrescos=cos((xres.').*freq);
xy=1*(pi/2).*sign(xrescos.*yrescos); %1 for maximum diffraction from 0 order
pattern_scatter=xy;
slm.show(pattern_scatter); 

%Get slm window data to update
slm.image_handle.Parent.CLim=[-pi,pi];
Cmap=slm.lookupTable.value;
Cmap=double(Cmap)./255;
slm.figure_handle.Colormap = Cmap;
slm.image_handle.CDataMapping = 'scaled';
slm.image_handle.CData = pattern_scatter;

%Power adjustment - see dynamic adjustmet in main loop
%(alternative/experimental)
load('slmLUT_191121.mat'); %%Correct this to account for ANGLES 
[Xq,Yq] = meshgrid(linspace(-max(rSpot),max(rSpot),dim)....
        ,linspace(-max(rSpot),max(rSpot),dim));
powerSurface =fittedmodel(Xq,Yq);
powerSurface = powerSurface./(max(max(powerSurface)));
powerSurface = smoothdata(powerSurface);
P0Max = min(min(powerSurface));

%% Setup multiple spots
assert(length(rSpot)<5 && length(rSpot)==length(thetaSpot),'Maximum of 4 diffracted spots allowed');
patternAll = {};
powerFactor = reducePowerBy + ((1-P0Max));
[~,ditheredPattern] = modulateAmplitude(powerFactor,freq,sz);
for nn=1:length(rSpot)
  %normalise spot power
  [xSpot,ySpot] = pol2cart(thetaSpot(nn),rSpot(nn));
  spotPower(nn) = fittedmodel(xSpot,ySpot); %power correction at spot point
  scale_value = dist_to_scale(rSpot(nn));
  pat = otslm.simple.linear(sz,scale_value,'angle',thetaSpot(nn)+offset);
  patternAll{nn} = pat;
end
pattern = otslm.tools.combine(patternAll,'method','rsuper','weights',spotPower);
pattern = otslm.tools.finalize(pattern);
slm.show(pattern+zeroCorrection)
%% Get 3d array from pattern
%Circle (cylinder), triangle, rect

%NOTE: include 'image_name','xxx' ...  if using bitmap generator
% [shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,...
%     'xdim',xlength,'ydim',ylength,'shape',structure,'image',im_name);
[shape_bitmap] = shape_map_2d_parse(dim,'coords',coords,...
       'xdim',xlength,'ydim',ylength,'shape','triangle');

[shape_bitmap_3d] = create_3d_structure(shape_bitmap,height_struct,res_struct);
%Display 3d bitmap to be printed
[~] = generate_3d_figure(shape_bitmap_3d);

% [shape_bitmap_3d]=generate_3d_sphere_grid(xlength,dim,0.5e-6);%,[-1e-6,0]);%(xlength,dim,900e-9); %Spheres

%% 
if ceil(height_struct/res_struct) >= 2
    sz_bitmap=size(shape_bitmap_3d);
    zsteps=sz_bitmap(3);
elseif ceil(height_struct/res_struct) == 1
    zsteps=1;
end

exposedPoints = length(find(squeeze(shape_bitmap_3d(:,:,:))==1));
timeEst = (dwelTime*exposedPoints)/60;
timeLayer = timeEst/layers;

disp(['Time estimate for 1 layer: ',num2str(timeLayer), ' Min.']);
disp(['Time estimate for structure: ',num2str(timeEst), ' Min.']);

disp('press any key to start...')
pause;
%% Start priting pattern
k=0;
firstRun = 1;
xgrid=linspace(coords(1),coords(1)+xlength,dim);
ygrid=linspace(coords(2),coords(2)+ylength,dim);
xgrid=linspace(coords(1),coords(1)+xlength,dim);
zgrid=linspace(0, res_struct*layers,layers);
counter =0;

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
end

initialVoltage = inputSingleScan(d);
for ii=1:layers
  z = zgrid(ii);
  for jj=1:dim
    y = ygrid(jj);
    for kk=1:dim
      
      x=xgrid(kk);
      [xv,yv,zv] = dist_to_volt_3d(x,y,z);
      pos = [xv,yv,zv];      
      
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
slm.image_handle.CData = pattern_scatter; 
disp(['Program Completed.']);

outputSingleScan(d,[0,0,0,]);

