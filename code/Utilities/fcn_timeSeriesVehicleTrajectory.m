%%%%%%%%%%%%%% Function fcn_timeSeriesVehicleTrajectory %%%%%%%%%%%%%%%%%%%
% Purpose:
%   fcn_timeSeriesVehicleTrajectory creates a timeseries object of 
%   swervedVehicleTrajectory
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   [timeSeriesVehicleTrajectory, sim_duration_from_aimsun_time] = ...
%   fcn_timeSeriesVehicleTrajectory(swervedVehicleTrajectory, fields, timeLookAhead)
% 
% INPUTS:
%   swervedVehicleTrajectory: queryVehicleTrajectory after adding sine
%   distortion
%   timeLookAhead: variable to decide how far to look ahead in time
%   fields: structure mapping column numbers of swervedVehicleTrajectory to
%   column names
% 
% OUTPUTS:
%   timeSeriesVehicleTrajectory: timeseries object
% 
% Author: Satya Prasad
% Create Date: 2020-05-29
% Revision history:
% 2020/05/30: Added the variable that indicates look ahead in time
% 2020/06/02: Added a structure containing column numbers as an input
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeSeriesVehicleTrajectory, sim_duration_from_aimsun_time] = ...
   fcn_timeSeriesVehicleTrajectory(swervedVehicleTrajectory, fields, timeLookAhead)
%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are there right number of inputs?
switch nargin
    case 3
    
    case 2
        timeLookAhead = 0;
    otherwise
        error('Incorrect number of input arguments');
end
% check the validity of timeLookAhead
if timeLookAhead < 0
    error('timeLookAhead must be non-negative number');
end
%% Convert swervedVehicleTrajectory into a timeseries object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_aimsun_time = max(swervedVehicleTrajectory(:,fields.aimsunTime));
min_aimsun_time = min(swervedVehicleTrajectory(:,fields.aimsunTime));
% duration of the simulation
sim_duration_from_aimsun_time = max_aimsun_time-min_aimsun_time;

% starting index of the trajectory
start_of_trajectory = round(10*timeLookAhead)+1;
% ending index of the trajectory
end_of_trajectory   = size(swervedVehicleTrajectory,1);
% length of the time object
length_of_time = end_of_trajectory-start_of_trajectory;

% time for timeseries object
timevals = 0.1*(0:length_of_time);
% create a time series object
timeSeriesVehicleTrajectory = ...
    timeseries(swervedVehicleTrajectory(start_of_trajectory:end_of_trajectory,:), timevals);

end