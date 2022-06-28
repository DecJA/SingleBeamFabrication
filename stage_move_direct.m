function [Vend] = stage_move_direct(Vstart,Vfinal)

%Find and assign DAQ board
%Locate DAQ board and channels
board = daq.getDevices;
inputs = board(1);

%Setup data recording session
sesh = daq.createSession('ni');

Vx0 = Vstart(1);
outputSignal = Vfinal(1);
%X-stage channel
[chx,idx] = addAnalogOutputChannel(sesh,'Dev2','ao0','Voltage');

disp('Starting Stage Movement');

queueOutputData(sesh,outputSignal);
startForeground(sesh);

%Update current postion
Vend = Vfinal;
disp('Movement Complete!');  


end


  