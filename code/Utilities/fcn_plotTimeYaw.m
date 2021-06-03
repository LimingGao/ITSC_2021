%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotTimeYaw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the time and yaw
%
% Input Variables:
%       tout: simualtion time
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
% fcn_plotTimeYaw(P_CG)
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
% Created Date:       2020-09-08
% Revisions:
%         
%
% To do list:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn_plotTimeYaw(tout,queryVehicleTrajectory,inputTrajectory,Desired_trajectory,P_CG)

h_fig = figure(753782);
set(h_fig,'Name','time and yaw');
clf;
hold on
plot(queryVehicleTrajectory(:,19)-queryVehicleTrajectory(1,19) ,queryVehicleTrajectory(:,15),'b-+','MarkerSize',6); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(inputTrajectory(:,19)-inputTrajectory(1,19),90-inputTrajectory(:,15)*180/pi,'c-o','MarkerSize',5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(tout,90- Desired_trajectory(:,7)*180/pi,'m-','MarkerSize',10); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(tout,90-P_CG(:,4)*180/pi,'r.','MarkerSize',10); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)

grid on 
legend('Reference Path from Aimsun','Input Trajectory to ST Trajectory query block','Desired Path yaw ','Simulated Path','Location','best');
title('Vehicle Yaw ');
xlabel('time [s]')
ylabel('yaw [degree]')
%xlim([min(station_cg)-10 1000+10])
box on
end

