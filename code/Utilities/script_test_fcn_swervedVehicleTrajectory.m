% fcn_swervedVehicleTrajectory adds sine wave distortion to 
% queryVehicleTrajectory iff flag_swerve is true
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   swervedVehicleTrajectory = fcn_swervedVehicleTrajectory(queryVehicleTrajectory,...
%                      amplitude, time_period, refLLA, fields, flag_swerve)
% 
% INPUTS:
%   queryVehicleTrajectory: Nx20 matrix where N is atleast 1.
% 
%   amplitude: It is either a row or column vector.
% 
%   time_period: It is either a row or column vector. amplitude and
%   time_period must have same length.
% 
%   refLLA: reference LLA for ENU to LLA transformation
% 
%   fields: structure mapping column numbers of queryVehicleTrajectory to column names
% 
%   flag_swerve: 'true' indicates add distortion. 'false' indicates no
%   distortion
% 
% OUTPUTS:
%   swervedVehicleTrajectory: Contains swerved trajectory
% 
% Author: Satya Prasad
% Create Date: 2020-05-28
% Revision history:
% 2020/06/02: Updated the script to test the modified
% fcn_swervedVehicleTrajectory
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
%% Only one sine wave
% swerve_amplitude = 1.7; swerve_timePeriod = 60; 
% flag.swerve = true;
% 
% swervedVehicleTrajectory = fcn_swervedVehicleTrajectory(queryVehicleTrajectory,...
%                             swerve_amplitude, swerve_timePeriod, refLLA, ...
%                             fieldsPreSimulationVehicleTrajectory, flag.swerve);

%% Two sine waves
swerve_amplitude = [1.7, 0.3]; swerve_timePeriod = [60, 5]; 
flag.swerve = true;

swervedVehicleTrajectory = fcn_swervedVehicleTrajectory(queryVehicleTrajectory,...
                            swerve_amplitude, swerve_timePeriod, refLLA, ...
                            fieldsPreSimulationVehicleTrajectory, flag.swerve);

%% Plots
reference_east  = queryVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.east);
reference_north = queryVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.north);

swerved_east  = swervedVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.east);
swerved_north = swervedVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.north);

figure(12345);
clf
plot(reference_east, reference_north, 'b*');
hold on;
plot(swerved_east, swerved_north, 'gd');
title('Vehicle Trajectory');
xlabel('east(m)'); ylabel('north(m)');
legend('Reference Trajectory', 'Swerved Trajectory');