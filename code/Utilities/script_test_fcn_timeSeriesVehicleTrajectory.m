% fcn_timeSeriesVehicleTrajectory creates a timeseries object of 
% swervedVehicleTrajectory
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
% 2020/06/02: Updated the script to test the modified
% fcn_timeSeriesVehicleTrajectory
% ======== to do list ============
% 1.   
%% Prepare the workspace
clear
clc

addpath('./Utilities/');
addpath('./DataFiles/');

load('vehicle_trajectory_2020_05_17.mat')
queryVehicleTrajectory = vehicle_trajectory.signals.values;
% reference LLA for ENU to LLA transformation
refLLA = [-48.876667, -123.393333, 0];
% field/column names of vehicle trajectory array before being used in the
% simulation
fieldsPreSimulationVehicleTrajectory.lat        = 9;
fieldsPreSimulationVehicleTrajectory.lon        = 10;
fieldsPreSimulationVehicleTrajectory.alt        = 11;
fieldsPreSimulationVehicleTrajectory.east       = 12;
fieldsPreSimulationVehicleTrajectory.north      = 13;
fieldsPreSimulationVehicleTrajectory.up         = 14;
fieldsPreSimulationVehicleTrajectory.yaw        = 15;
fieldsPreSimulationVehicleTrajectory.speed      = 16;
fieldsPreSimulationVehicleTrajectory.station    = 18;
fieldsPreSimulationVehicleTrajectory.aimsunTime = 19;
% add swerving to the queryVehicleTrajectory
swerve_amplitude = [1.7, 0.3]; swerve_timePeriod = [60, 5]; 
flag.swerve = true;

swervedVehicleTrajectory = fcn_swervedVehicleTrajectory(queryVehicleTrajectory,...
                            swerve_amplitude, swerve_timePeriod, refLLA, ...
                            fieldsPreSimulationVehicleTrajectory, flag.swerve);

%% create time series object for vehicle trajectory
timeLookAhead = 1;
[timeSeriesVehicleTrajectory, sim_duration_from_aimsun_time] = ...
        fcn_timeSeriesVehicleTrajectory(swervedVehicleTrajectory, ...
        fieldsPreSimulationVehicleTrajectory, timeLookAhead); 

%% Run the simulation
% simulink model takes in a time series object and outputs the same
% trajectory interpolated based on the simulation step size
sim('simulink_test_fcn_timeSeriesVehicleTrajectory.slx', sim_duration_from_aimsun_time);

%% Plot for both input and output trajectory
swerved_east  = swervedVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.east);
swerved_north = swervedVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.north);

input_east  = timeSeriesVehicleTrajectory.Data(:,fieldsPreSimulationVehicleTrajectory.east);
input_north = timeSeriesVehicleTrajectory.Data(:,fieldsPreSimulationVehicleTrajectory.north);

output_east  = output_trajectory(:,4);
output_north = output_trajectory(:,5);

figure(12345);
clf;
plot(swerved_east,swerved_north,'b.');
hold all
plot(input_east,input_north,'rd','MarkerSize',12);
plot(output_east,output_north,'go');
title('Vehicle Trajectory');
legend('Swerved Trajectory','Input Trajectory','Output Trajectory');
xlabel('East (m)'); ylabel('North (m)');