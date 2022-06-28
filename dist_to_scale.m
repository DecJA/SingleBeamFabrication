function [scale_input,params] = dist_to_scale(desired_x)
%For a specified disance from oth order, calculates input for otslm
%function

%Power fit parameters calculated through R:\OMG\declan\TableD\slmShift
%CCam = 1.1712e-07; %1.2854e-07; %10/05/2021
% CCam = 1.0624e-07; %16/06/2021
% CCam = 1.075e-07; %30/09
% CCam = 1.108e-07; %12/11/21
% CCam = 1.424e-7; % 10/02
% CCam = 2.1322e-7; % 09052022
CCam = 1.579e-7; %20/06/22
%% Paramaters measured 30/09/2021
% A = 6.338e-5;
% B = -1; %-0.98 fit function
%% Measured 21/01/22
% A = 2.04e-4;%1.552e-4;
% B = -1;
% C = 0;
% scale_input = A*(desired_x^B) + C; %1.0 R2 value

% A =  1.9590e-04; %1.552e-4;
% B = -1;
% C = 0;

% A = 9.88e-5;
% B = -1;
% C = 0;

%Measured 20/06/22 - power fit from scale(X) vs colDist (col-min(col))*CCam
A = 1.492e-4;
B = -1;
C = 0;
scale_input = A*(desired_x.^B) + C; %1.0 R2 value

end