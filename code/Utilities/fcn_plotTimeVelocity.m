%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotTimeVelocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the time and velocity
%
% Input Variables:
%       tout: simualtion time
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
% fcn_plotTimeVelocity(tout,)
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

function fcn_plotTimeVelocity(tout,command_velocity,inputTrajectory,Desired_trajectory,Vx,Vy,safe_velocity)
            
h_fig = figure(726982);
width=600;%
height=450;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
set(h_fig,'Name','time and velocity');
clf;
hold on
safe_velocity(safe_velocity>=36) = 36;
%plot(inputTrajectory(:,19)-inputTrajectory(1,19),inputTrajectory(:,16),'c-o','MarkerSize',5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(tout,Desired_trajectory(:,8),'c--','MarkerSize',10,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(tout,sqrt(Vx.^2+Vy.^2),'m.','MarkerSize',10,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(tout,safe_velocity,'r-','MarkerSize',10,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)
plot(tout,command_velocity,'b-','MarkerSize',6,'LineWidth',1.5); % convert the yaw definition from atan2 to IVSG(%Yaw,north is zero, clockwise is positive direction, range 0-360 degrees)

grid on 
legend('Desired Tracking Speed ','Simulated Speed','Previewed Allowable Speed','Command Speed','Location','best');
% title('Vehicle Speed ');
xlabel('Time [s]')
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

