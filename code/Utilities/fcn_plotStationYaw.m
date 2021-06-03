%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationYaw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the tation of cg postion and yaw
%
% Input Variables:
%       queryVehicleTrajectory:vehicle path from AIMSUN
%       inputTrajectory: input vehicle path to ST Trajectory Query block
%       Desired_trajectory: output vehicle path from  ST Trajectory Query block
%       % P_CG_ENH: output vehicle path from vehicle model
%       P_CG: the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
% 
% 
% Returned Results:
%
% Example:
% fcn_plotStationYaw(P_CG)
% 
% Processing Flow:
%  IVSG yaw defintion: north is zero, clockwise is positive direction, range 0-360 degrees
%
% Restrictions/Notes:
%
% The following functions are called:
%      none
%
% Author:             Liming Gao
% Created Date:       2020-05-28
% Revisions:
%           2020-05-29: 
%
% To do list:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn_plotStationYaw(queryVehicleTrajectory,inputTrajectory,Desired_trajectory,P_CG)

% station_aimsun = [0; cumsum(sqrt(sum(diff(queryVehicleTrajectory(:,12:13)).^2,2)))];
station_inputTrajectory = [0; cumsum(sqrt(sum(diff(inputTrajectory(:,12:13)).^2,2)))];
station_DesiredTrajectory = [0; cumsum(sqrt(sum(diff(Desired_trajectory(:,4:5)).^2,2)))];
station_cg = [0; cumsum(sqrt(sum(diff(P_CG(:,1:2)).^2,2)))]; % calculate station

h_fig = figure(75389);
set(h_fig,'Name','station of cg postion and yaw');
clf;
hold on
plot(queryVehicleTrajectory(:,18),queryVehicleTrajectory(:,15),'b-+','MarkerSize',6); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(inputTrajectory(:,18),90-inputTrajectory(:,15)*180/pi,'c-o','MarkerSize',5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(station_cg,90- Desired_trajectory(:,7)*180/pi,'m-','MarkerSize',10); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(station_cg,90-P_CG(:,4)*180/pi,'r.','MarkerSize',10); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)

grid on 
legend('Reference Path from Aimsun','Input Trajectory to ST Trajectory query block','Desired Path yaw ','Simulated Path','Location','best');
title('Vehicle Yaw ');
xlabel('station [m]')
ylabel('yaw [degree]')
%xlim([min(station_cg)-10 1000+10])
box on
end

