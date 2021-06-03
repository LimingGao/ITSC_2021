%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationWheelTorque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      plot the tation of cg postion and wheel driving torque
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      wheel_torques = driving torque of four wheels, format: N*4 matrix
%      varargin{1} =  driving torque of four wheels without friction preview, format: N*4 matrix

% Returned Results:
%
% Example:
%  fcn_plotStationWheelTorque(P_CG,wheel_torques,wheel_torques_without_preview)
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
% Created Date:       2020-09-18
% Revisions:
%          
%
% To do list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn_plotStationWheelTorque(P_CG,wheel_torques,varargin)

h_fig = figure(19432);
set(h_fig,'Name','station and wheels torque');
width=600;%
height=450;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on
station_cg = [0; cumsum(sqrt(sum(diff(P_CG(:,1:2)).^2,2)))]; % calculate station
if nargin ==2
    plot(station_cg,wheel_torques(:,4),'k-','MarkerSize',6,'LineWidth',2);
    legend_text ={'rear right wheel driving torque'};
    
elseif nargin ==3
    plot(station_cg,wheel_torques(:,4),'k-','MarkerSize',6,'LineWidth',2);
    wheel_torques_without_preview = varargin{1};
    plot(station_cg,wheel_torques_without_preview(:,4),'b--','MarkerSize',6,'LineWidth',2.3);
    legend_text ={'with friction preview','without friction preview'};
    
end

grid on
legend(legend_text,'Location','Best');
% title('Wheels Torque v.s. Station ');
xlabel('Station [m]','fontsize',14)
ylabel('Wheels Torque [Nm]','fontsize',14)
ylim([-250 250])
%xlim([0 1050])
axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
axes_handle.GridLineStyle = '-.';
axes_handle.GridColor = 'k';
axes_handle.GridAlpha = 0.2;
%             ax.YMinorGrid = 'on';
%             ax.MinorGridAlpha = 0.5
box on


end

