function [board_check,slm_check] = system_run(dim,coords,xlength,ylength,zlength,structure,layers)

offset = -(3/2)*pi; %Offset from SLM-Camera angle


%% Check if SLM is setup yet
slm_check = exist('sz','var');
if slm_check==0  
    addpath('R:\declan\TableD');
    addpath('R:\declan\TableD\Printing');
    cd('R:\declan\TableD');
    run('setup_slm.m');
    cd('R:\declan\TableD\Printing\CoherentProgram');
    slm_check = 1;
end

%% Check if DAQ board is connected
board_check = exist('board_check','var');
if board_check==0
    board = daq.getDevices;
    inputs = board(1);

    %Get current voltage output from channels
    sesh = daq.createSession('ni');
    sesh.DurationInSeconds = 1.0; %Record one second of data
    addAnalogInputChannel(sesh,'Dev2',0,'Voltage');
    board_check = 1;
end

%% Spectify grid area and print resoltion

%Axial steps and voltages
res_struct = zlength;
height_struct = res_struct*layers;  
step_voltage = dist_to_volt(res_struct);

%display scattering pattern on SLM to avoid polymerisation
pattern_scatter = otslm.simple.checkerboard(sz);
pattern_scatter = otslm.tools.finalize(pattern_scatter);
slm.show(pattern_scatter); 

[shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,'xdim',xlength,'ydim',ylength,'shape',structure);
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
    
    %Get current axial voltage
    [data,time] = sesh.startForeground;
    V0 = mean(data);
    if k==0
        V0N = V0;
        k = k+1;
    else
        %Send new voltage to device
        [Vx_act] = stage_move_direct(V0,V0-step_voltage);
        %Confirm voltage target    
        [dataN,time] = sesh.startForeground;
        V0N = mean(dataN);
        %increment voltage step
        k = k+1;
    end
    
    % Lateral SLM step
    
    for ii=1:dim
        x=xgrid(ii);
        for jj=1:dim
            y=ygrid(jj);
            
            %Check if print required from bitmap
            bit_point = bitmap_plane(ii,jj);
            if bit_point == 1
                
                [theta,r] = cart2pol(x,y);
                scale_value = dist_to_scale(r);
                pattern = otslm.simple.linear(sz,scale_value,'angle',theta+offset);
                pattern = otslm.tools.finalize(pattern);
                slm.show(pattern);
                %pause(1/60);
                
            elseif bit_point == 0
                %slm.show(pattern_scatter);
                continue
            end
            
        end
    end
                
    disp(['Moving to step: ',num2str(k),' of', num2str(zsteps-1)]);
    
end
slm.show(pattern_scatter);
disp(['Program Completed.']);

end