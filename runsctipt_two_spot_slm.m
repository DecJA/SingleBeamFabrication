
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
startAxialVoltage = 5;
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
xlength = 8e-06;
ylength = 8e-06;

%Axial steps and voltages
res_struct = 200e-09;
layers =5;

structure = 'image';%Botton left grid cornder locdddaation
im_name = 'barrier_array2';
coords= [2e-6,0e-06];
dim = 50;

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

%Sample power adjustment
reducePowerBy = 0; %0-no reduction - 1 - 100% power reduceion (works up to 0.90)
if reducePowerBy ~=0
    [scaleFactor,ditheredPattern] = modulateAmplitude(reducePowerBy,128,sz); %50% power modulation
    slm.image_handle.CData = ditheredPattern;
else
    ditheredPattern = zeros(sz);
end

%% Get 3d array from pattern
%Circle (cylinder), triangle, rect

%NOTE: include 'image_name','xxx' ...  if using bitmap generator
[shape_bitmap] = shape_map_2d_parse(dim,'coords',coords,...
       'xdim',xlength,'ydim',ylength,'shape','triangle');
%[shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,...
  %   'xdim',xlength,'ydim',ylength,'shape',structure,'image',im_name);
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
k=0;
xgrid=linspace(coords(1),coords(1)+xlength,dim);
ygrid=linspace(coords(2),coords(2)+ylength,dim);

while k<zsteps
    %Get slice from structure
    bitmap_plane = shape_bitmap_3d(:,:,k+1);
    %Get current axial voltage and move if needed
    currentVoltage=inputSingleScan(d);
    if k~=0
        disp(['Moving to step: ',num2str(k),' of', num2str(zsteps)]);
        outputSingleScan(d,currentVoltage-step_voltage);
    end
    
    % Lateral SLM step -- linear indexing
    lin_idx = find(bitmap_plane>0);
    lin_idx1 = lin_idx(1:end/2);
    lin_idx2 = flip(lin_idx(end/2+1:end));
    
    assert(length(lin_idx1) == length(lin_idx2),...
              'WARNING!: Linear Index Mismatch!');
    
    for j=1:numel(lin_idx1)
        
        [row1,col1]=ind2sub(size(bitmap_plane),lin_idx1(j));
        [row2,col2]=ind2sub(size(bitmap_plane),lin_idx2(j));
        
         x1=xgrid(row1); y1=ygrid(col1); 
         x2=xgrid(row2); y2=ygrid(col2); 

         [theta1,r1] = cart2pol(x1,-y1);
         [theta2,r2] = cart2pol(x2,-y2);
         scale_value1 = dist_to_scale(r1);
         scale_value2 = dist_to_scale(r2);
         
         sp1 = otslm.simple.linear(sz,scale_value1,'angle',theta1+offset);
         sp2 = otslm.simple.linear(sz,scale_value2,'angle',theta2+offset);
         pattern = otslm.tools.combine({sp1,sp2},'method','super');
         
         %To Do: include distance-power modulation
         pattern = otslm.tools.finalize(pattern);
         allPatterns = pattern + zeroCorrection + ditheredPattern;
         slm.image_handle.CData = allPatterns;
         pause(0.1);
    end
    slm.image_handle.CData = pattern_scatter; 
    k=k+1;
    pause(0.5);   
end

disp(['Program Completed.']);


