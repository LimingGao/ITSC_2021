%%%%%%%%%%%%%%%%%%% Function fcn_plotTrajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcn_plotTrajectory plots vehicle path from AIMSUN, ST Trajectory Query
% block and Vehicle Dynamics simulation
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
%
% Format:
%   fcn_plotTrajectory(eastTraj, northTraj)
% 
% INPUTS:
%   queryVehicleTrajectory:vehicle path from AIMSUN
%   inputTrajectory: input vehicle path to ST Trajectory Query block
%   Desired_trajectory: output vehicle path from  ST Trajectory Query block
%   P_CG_ENH: output vehicle path from vehicle model
%   preview_friction: the coordinate of preview point
% 
% OUTPUTS:
%   plot of trajectories and yaw angle 
% 
% Author: Satya Prasad, Liming
% Create Date: 2020-05-29
% Revision history:
% 2020/06/02: Added the visualization plot for vehicle trajectory
% 2020/06/04: Modified inputs from arrays to struct
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  fcn_plotTrajectory(queryVehicleTrajectory,inputTrajectory,Desired_trajectory,P_CG_ENH,preview_friction)
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
if nargin ~= 5 && nargin ~= 1
    error('Incorrect number of input arguments');
end

%% plot vehicle trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% east and north coordinates on vehicle trajectory in the database/from AIMSUN
trajEast.reference  = queryVehicleTrajectory(:,12);
trajNorth.reference = queryVehicleTrajectory(:,13);

if nargin == 5
    % input of ST Trajectory DB Query
    trajEast.inputTrajectory = inputTrajectory(:,12);
    trajNorth.inputTrajectory = inputTrajectory(:,13);
    
    % east and north coordinates on vehicle trajectory from the ST Trajectory Query block
    trajEast.desired  = Desired_trajectory(:,4);
    trajNorth.desired = Desired_trajectory(:,5);
    
    % east and north coordinates on vehicle trajectory, ouput of the simulation
    trajEast.simulated  = P_CG_ENH(:,1);
    trajNorth.simulated = P_CG_ENH(:,2);
    
    % eastTraj: east coordinates on trajectory
    % northTraj: north coordinates on trajectory
    
    h_fig = figure(512345);
    set(h_fig,'Name','compare vehicle trajectories');
    clf;
    plot(trajEast.reference,trajNorth.reference,'r');
    hold all;
    plot(trajEast.inputTrajectory,trajNorth.inputTrajectory,'c:','LineWidth',2);
    plot(trajEast.desired,trajNorth.desired,'b*');
    plot(trajEast.simulated,trajNorth.simulated,'gd');
    plot(preview_friction(:,1),preview_friction(:,2),'m-s','MarkerSize',10);
    legend('Reference Path from Aimsun','Input Trajectory to ST Trajectory query block','Desired Path','Simulated Path','Preview Location','Location','best');
    title('Vehicle Trajectory');
    xlabel('East (m)'); ylabel('North (m)');
    grid on
elseif nargin == 1
    h_fig = figure(513345);
    set(h_fig,'Name','fcn_plotTrajectory:queryVehicleTrajectory');
    clf;
    plot(trajEast.reference,trajNorth.reference,'r');
    
    legend('Reference Path from Aimsun','Location','best');
    title('Vehicle Trajectory');
    xlabel('East (m)'); ylabel('North (m)');
    grid on
else
    
     error('Incorrect number of input arguments');
    
end
end