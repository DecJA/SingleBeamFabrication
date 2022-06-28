%% V2 updates
%Corrected circle shift function
%Make stage move in one 1step (still need to test)
%allow for 2d shape printing

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
  [chi,idi] = addAnalogInputChannel(d,'Dev2',[0:3],'Voltage');
  %removeChannels(d,1:3);
  [ch,id] = addAnalogOutputChannel(d,'Dev2',[0:3],'Voltage'); %Add 
  outputSingleScan(d,[0,0,5.5,1]);
  initialVoltage = inputSingleScan(d);
   board_check = 1;    
end

%% Spectify grid area and print resoltion
cur=pwd;
maps=strcat(cur,'\PrintMaps');
addpath(maps);

offset =-(3/2)*pi; %Offset from SLM-Camera angle
xlength = 8e-06;
ylength = 8e-06;

%Axial steps and voltages
res_struct = 200e-09;
layers =20;

structure = 'image';%Botton left grid cornder locdddaation
im_name = 'barrier_simple2';
coords= [0,0];
dim = 80;

voxel_spacing = (max(xlength)/dim)*1e9;
c=newline;
disp([c,'Voxel spacing is set to: ',num2str(voxel_spacing),' nm.']);

dwelTime = 0.05;
reducePowerBy = 0.5; %0-no reduction - 1 - 100% power reduceion (works up to 0.90)

%% Setup axial step and slm
height_struct = res_struct*layers;  
step_voltage = dist_to_volt(res_struct);

%Get current stage voltage and reposition
% currentVoltage=inputSingleScan(d);
% if currentVoltage(4) < 0.1 %Check shutter position
%   outputSingleScan(d,[0,0,5.5,1]); %Close Shutter
% end
% if round(currentVoltage(3),3)~=round(initialVoltage(3),3) %check axial height
%   outputSingleScan(d,initialVoltage);
% end
outputSingleScan(d,initialVoltage);
%%  display scattering pattern on SLM to avoid polymerisation
xres = linspace(-1,1, sz(1));
freq=128;
xrescos=cos(xres.*(freq));
yrescos=cos((xres.').*freq);
xy=1*(pi/2).*sign(xrescos.*yrescos); %1 for maximum diffraction from 0 order
pattern_scatter=xy;
slm.show(pattern_scatter); 

%Sample power adjustment
if reducePowerBy ~=0
    [~,ditheredPattern] = modulateAmplitude(reducePowerBy,128,sz); 
    [~,ditheredPatternMid] = modulateAmplitude(reducePowerBy+0.15,128,sz);
    slm.show(ditheredPattern);
else
    ditheredPattern = zeros(sz);
end

%% Get 3d array from pattern
%Circle (cylinder), triangle, rect
%NOTE: include 'image_name','xxx' ...  if using bitmap generator
% [shape_bitmap] = shape_map_2d_parse(dim,'coords',coords,...
%        'xdim',xlength,'ydim',ylength,'shape','triangle');
[shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,...
    'xdim',xlength,'ydim',ylength,'shape',structure,'image',im_name);
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

Ntotal = layers*dim*dim;

exposedPoints = length(find(squeeze(shape_bitmap_3d(:,:,1))==0));
timeLayer = (dwelTime*exposedPoints)/60;
timeEst = (timeLayer*layers);

disp(['Time estimate for 1 layer: ',num2str(timeLayer), ' Min.']);
disp(['Time estimate for structure: ',num2str(timeEst), ' Min.']);

%Get slm window data to update
slm.image_handle.Parent.CLim=[-pi,pi];
Cmap=slm.lookupTable.value;
Cmap=double(Cmap)./255;

slm.figure_handle.Colormap = Cmap;
slm.image_handle.CDataMapping = 'scaled';
slm.image_handle.CData = pattern_scatter;

%%
for ii=1:layers
  z = zgrid(ii);
  if ii==1
    patternCombined = ditheredPattern+zeroCorrection+testHex;
    slm.image_handle.CData = patternCombined;
    %outputSingleScan(d,[currentVoltage(1),currentVoltage(2),currentVoltage(3),0]);
  else
    patternCombined = ditheredPatternMid+zeroCorrection;
    slm.image_handle.CData = patternCombined;  
    %outputSingleScan(d,[pos(1),pos(2),pos(3),0]);
  end
  for jj=1:dim
    y=ygrid(jj);
    for kk=1:dim
      
      x=xgrid(kk);
      [xv,yv,zv] = dist_to_volt_3d(x,y,z);
      pos = initialVoltage+[xv,yv,-zv,0];      
      
      if shape_bitmap_3d(kk,jj,ii)==1
        outputSingleScan(d,pos);
        pause(dwelTime);
      end 
    end
  end
  %outputSingleScan(d,[pos(1),pos(2),pos(3),1]); %Close Shutter
  slm.image_handle.CData = pattern_scatter;
  disp(['Moving to step: ',num2str(ii),' of ', num2str(zsteps)]);
end
outputSingleScan(d,initialVoltage); %Close Shutter
disp(['Program Completed.']);

