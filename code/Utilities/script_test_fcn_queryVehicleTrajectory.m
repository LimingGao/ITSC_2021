% fcn_queryVehicleTrajectory queries a vehicle's trajectory uniquely 
% identified by (vehicle_id, global_time)
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, global_time)
% 
% INPUTS:
%   vehicle_id: integer of size one. ID of the vehicle
%   global_time: intger of size one. Uniquely identifies a vehicle's 
%   trajectory in the database
% 
% OUTPUTS:
%   queryVehicleTrajectory: Contains all the attributes defined by
%   trajectory_attributes sorted in the order of aimsun_time
% 
% Author: Satya Prasad
% Create Date: 2020-05-27
% Revision history:
% 1. 
% ======== to do list ============
% 1. 
%% Prepare the workspace
clear
clc

addpath('./Utilities/');
addpath('./DataFiles/');

%% Show error with 'vehicle_id'
% vehicle_id = -1; global_time = 1.589072708000000e+09;
% queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, global_time);

%% Show error with 'global_time'
% vehicle_id = 26; global_time = -1;
% queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, global_time);

%% Query a valid trajectory
vehicle_id = 26; global_time = 1.589072708000000e+09;
queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, global_time);