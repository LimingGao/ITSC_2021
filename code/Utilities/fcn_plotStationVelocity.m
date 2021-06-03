%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationVelocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the tation of cg postion and velocity
%
% Input Variables:
%       command_velocity:velocity command to controller
%       inputTrajectory: input vehicle path to ST Trajectory Query block
%       Desired_trajectory: output vehicle path from  ST Trajectory Query block
%       % P_CG_ENH: output vehicle path from vehicle model
%       P_CG: the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%       safe_velocity: maximum allowable velocity
% 
% Returned Results:
%
% Example:
% fcn_plotStationVelocity(P_CG)
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
% To do list:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   fcn_plotStationVelocity(command_velocity,inputTrajectory,Desired_trajectory,Vx,Vy,ahead_safe_velocity,P_CG)
            
% station_aimsun = [0; cumsum(sqrt(sum(diff(queryVehicleTrajectory(:,12:13)).^2,2)))];
station_inputTrajectory = [0; cumsum(sqrt(sum(diff(inputTrajectory(:,12:13)).^2,2)))];
station_DesiredTrajectory = [0; cumsum(sqrt(sum(diff(Desired_trajectory(:,4:5)).^2,2)))];
station_cg = [0; cumsum(sqrt(sum(diff(P_CG(:,1:2)).^2,2)))]; % calculate station

h_fig = figure(7282);
width=600;%
height=450;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
set(h_fig,'Name','station and velocity');
clf;
hold on
ahead_safe_velocity(ahead_safe_velocity>=36) = 36;
%plot(inputTrajectory(:,19)-inputTrajectory(1,19),inputTrajectory(:,16),'c-o','MarkerSize',5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(station_cg,Desired_trajectory(:,8),'c--','MarkerSize',10,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(station_cg,sqrt(Vx.^2+Vy.^2),'m.','MarkerSize',10,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(station_cg,ahead_safe_velocity,'r-','MarkerSize',10,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(station_cg,command_velocity,'b-','MarkerSize',6,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)

grid on 
legend('Desired Tracking Speed ','Simulated Speed','Previewed Allowable Speed','Command Speed','Location','best');
% title('Vehicle Speed ');
xlabel('Station [m]')
ylabel('Speed [m/s]')
ylim([0 38])
%xlim([min(station_cg)-10 1000+10])
axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
axes_handle.GridLineStyle = '-.';
axes_handle.GridColor = 'k';
axes_handle.GridAlpha = 0.2;
%             ax.YMinorGrid = 'on';
%             ax.MinorGridAlpha = 0.5
box on
end

