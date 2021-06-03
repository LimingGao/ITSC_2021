function fcn_plotVehicleOnTrajectory(input_trajectory,E,N,PSI,delta,...
                                    Vehicle_Geometry,tout)
% fcn_plotVehicleOnTrajectory plots vehicle on the trajectory
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
%
% Format:
%   fcn_plotVehicleOnTrajectory(input_trajectory,E,N,PSI,delta,...
%                                     Vehicle_Geometry,tout)
% 
% INPUTS:
%   input_trajectory: vehicle trajectory input to the simulation
%   E,N,PSI: pose output of the vehicle dynamics simulation
%   delta: steering angle
%   Vehicle_Geometry: geometry of the vehicle
%   tout: simulation time
% 
% OUTPUTS:
%   plot of vehicle on trajectory
% 
% Author: Satya Prasad
% Create Date: 2020-05-29
% Revision history:
% 2020/02/06: Added the visualization plot for vehicle trajectory
% ======== to do list ============
% 1. 
% east and north coordinate input to the simulation
input_east  = input_trajectory(:,12);
input_north = input_trajectory(:,13);

hOH = figure(512346);
hold off;
plot(input_east,input_north,'k-.');
hold on;
vehicle.a = Vehicle_Geometry(1);
vehicle.b = Vehicle_Geometry(2);
vehicle.d = Vehicle_Geometry(3);
vehicle.Re = 0.32;
OverheadPlot(hOH,tout,E,N,PSI,delta,linspace(tout(1),tout(end),40),vehicle)
legend('Input Trajectory to ST Trajectory query block','Simulated Path','Location','best');

end