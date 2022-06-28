
slm_check = exist('sz','var');
if slm_check==0  
    addpath('R:\declan\TableD');
    addpath('R:\declan\TableD\Printing');
    cd('R:\declan\TableD');
    run('setup_slm.m');
    cd('R:\declan\TableD\Printing\CoherentProgram');
    slm_check = 1;
% %     load('zeroCorrectionSmooth.mat');
%     zeroCorrection = -zeroCorrectionSmooth';
    load('zeroCorrection.mat');
end

%% Check if DAQ board is connected
board_check = exist('board_check','var');
startAxialVoltage = 0;
if board_check==0
  d = daq.createSession('ni');
  %Make sure no voltage is output
  [chi,idi] = addAnalogInputChannel(d,'Dev2',2,'Voltage'); %Only z-axis
  [ch,id] = addAnalogOutputChannel(d,'Dev2',2,'Voltage'); %Add only z-axis
  outputSingleScan(d,startAxialVoltage); %Start with large positive voltage for axial dim.
  initialVoltage = inputSingleScan(d);
  board_check = 1;  
end

%% Spectify grid area and print resoltion
cur=pwd;
maps=strcat(cur,'\PrintMaps');
addpath(maps);

offset =-(3/2)*pi; %Offset from SLM-Camera angle
xlength = 14e-06;
ylength = 14e-06;

%Axial steps and voltages
layers = 3;
res_struct = 200e-09;
dPz=(linspace(0,.2.^1,layers+1));
dPz(1)=-.2;
dPz(2)=-.1;
dPz(7:end)=dPz(7);
dwelTime =0.04;
reducePowerBy = 0.5; %0-no reduction - 1 - 100% power reduceion (works up to 0.90)


structure = 'image';%Botton left grid cornder locdddaation
im_name = 'flower_new2';
coords= [-7e-06,-7e-06];
dim =80;
voxel_spacing = (max(xlength)/dim)*1e9;
c=newline;
disp([c,'Voxel spacing is set to: ',num2str(voxel_spacing),' nm.']);
%% Setup axial step and slm
height_struct = res_struct*layers;  
step_voltage = dist_to_volt(res_struct);

maxAxialVoltage = startAxialVoltage*(res_struct*step_voltage);
if maxAxialVoltage >= (startAxialVoltage + 2)
  warning('WARNING: Total axial displacement exceedes travel limit.');
end

%Get current stage voltage and reposition
currentVoltage=inputSingleScan(d);
if round(currentVoltage,1)~=round(initialVoltage)
  outputSingleScan(d,initialVoltage);
end

%%  display scattering pattern on SLM to avoid polymerisation
%   calculate dithering pattern to scale power
xres = linspace(-1,1, sz(1));
freq=128;
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
P0Max = translation_decrease_multiSpot(max(max([xlength; ylength])));
%% Get 3d array from pattern
%Circle (cylinder), triangle, rect

%NOTE: include 'image_name','xxx' ...  if using bitmap generator
[shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,...
    'xdim',xlength,'ydim',ylength,'shape',structure,'image',im_name);
% [shape_bitmap] = shape_map_2d_parse(dim,'coords',coords,...
%        'xdim',xlength,'ydim',ylength,'shape','triangle');

[shape_bitmap_3d] = create_3d_structure(shape_bitmap,height_struct,res_struct);
%Display 3d bitmap to be printed
[~] = generate_3d_figure(shape_bitmap_3d);

% [shape_bitmap_3d]=generate_3d_sphere_grid(xlength,dim,1e-6); %Spheres

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
counter =0;
while k<zsteps
    %Get slice from structure
    bitmap_plane = shape_bitmap_3d(:,:,k+1);
    %Get current axial voltage and move if needed
    currentVoltage=inputSingleScan(d);
    if k~=0
        point_check = length(find(bitmap_plane > 0));
        if point_check>1
            disp(['Moving to step: ',num2str(k),' of ', num2str(zsteps)]);
            outputSingleScan(d,currentVoltage+step_voltage);
            pause(1.0);
        end
    end
    
    % Lateral SLM step -- linear indexing
    lin_idx = find(bitmap_plane>0);
    for j=1:numel(lin_idx)
      
        tic
        [row,col]=ind2sub(size(bitmap_plane),lin_idx(j));
         x=xgrid(row);
         y=ygrid(col); 

         [theta,r] = cart2pol(x,-y);
         scale_value = dist_to_scale(r);
         
         pattern = otslm.simple.linear(sz,scale_value,'angle',theta+offset);
         pattern = otslm.tools.finalize(pattern);
         
         tempScaling =primary_power_scaling(P0Max,r);
         powerFactor = (1-tempScaling)+reducePowerBy;
         
         if powerFactor<0
           powerFactor = 0;
         elseif powerFactor>0.9
           disp(['Requested power reduction is: ',num2str(powerFactor), ...
                  ' Maximum possible dithering is ~ 90%']);
           powerFactor = 0.9;
         end
         [~,ditheredPattern] = modulateAmplitude(powerFactor*(1+dPz(k+1)),128,sz);
         allPatterns = pattern + zeroCorrection + ditheredPattern;
         slm.image_handle.CData = allPatterns;
         t2=toc;
         if dwelTime>(t2)
          pause(dwelTime-t2);
         end
    end
   
    slm.image_handle.CData = pattern_scatter; 
    pause(0.5);   
    if firstRun == 1
      firstRun =0;
    else
     k=k+1;
    end
end

slm.image_handle.CData = pattern_scatter; 
disp(['Program Completed.']);
%Reset stage position
outputSingleScan(d,currentVoltage);

%%
addAnalogOutputChannel(d,'Dev2',3,'Voltage'); %shutter
timer=[0:1:50];
for tt=1:length(timer)
  send =rem(timer(tt),2);
  outputSingleScan(d,[currentVoltage,send]);
  pause(0.05);
end
removeChannel(d,3)

