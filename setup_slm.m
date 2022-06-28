% Setup the SLM for general use
% clear all; clc;
addpath('R:\MatlabLibraries');
addpath('R:\EiffelTower');

%% Load the lookup table
% This one was exported from LabView but we could generate our own with
% the OTSLM calibration functions or just guessing the values.

fname = 'lookup_table_785_supplied.txt';%'SLM_colormap_25032019.txt';
%fname = 'small_lookup_table.txt';
lookup_table = otslm.utils.LookupTable.load(fname, ...
  'channels', [2, 2, 0], 'phase', [], 'format', @uint16, ...
  'mask', [hex2dec('00ff'), hex2dec('ff00')], 'morder',  1:8);

%% Setup the prescaled option for ScreenDevice.show
% The ScreenDevice.show function in the SLM takes the modulo of pattern
% values, to do this it needs to know what scaling the pattern values have

% Choose one of these options
% value_range = '0 to 2pi';
value_range = '-pi to pi'; %used 2
% value_range = '0 to 1'; %previously used
switch value_range
  case '0 to 2pi'
    prescaled = true;
  case '-pi to pi'
    prescaled = true;
  case '0 to 1'
    prescaled = false;
  otherwise
    error('Invalid option selected');
end

%% Setup the SLM
% Finally, we can setup the SLM

% Specify the device size in pixels
sz = [510, 510];

% Specify the monitor number we are using
n = 1;

% The monitor representing the device is larger than the device
% So we get the monitor size and place the device image in one corner
set(0,'units','pixels');
scsz = get(0,'MonitorPositions');  % Monitor sizes
screen_size = scsz(n, :);

% And finally, we get to creating the SLM object
slm = otslm.utils.ScreenDevice(n, 'target_size', sz, ...
  'target_offset', [0, screen_size(4)-sz(2)+1], ...
  'lookup_table', lookup_table, ...
  'pattern_type', 'phase', 'prescaledPatterns', prescaled);


%Check framereate setting in ScreenDevice function
% slm = otslm.utils.ScreenDevice(n, 'target_size', sz, ...
%   'target_offset', [0, screen_size(4)-sz(2)+1], ...
%   'lookup_table', lookup_table, ...
%   'pattern_type', 'phase', 'prescaledPatterns', prescaled,...
%   'value_range',{[0:255]});
%%
pattern_grate = otslm.simple.linear(slm.size,70);
pattern = otslm.simple.lgmode(sz,3,0);
pattern3 = otslm.tools.finalize(pattern+pattern_grate);
slm.show(pattern3)