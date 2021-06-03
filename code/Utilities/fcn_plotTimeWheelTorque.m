%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotTimeWheelTorque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      plot the time and wheel driving torque
%
% Input Variables:
%      tout   = time step,format: N*1 matrix
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

function fcn_plotTimeWheelTorque(tout,wheel_torques,varargin)

h_fig = figure(19433);
set(h_fig,'Name','time and wheels torque');
width=600;%
height=450;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on

if nargin ==2
    plot(tout,wheel_torques(:,4),'k-','MarkerSize',6,'LineWidth',2);
    legend_text ={'rear right wheel driving torque'};
    
elseif nargin ==3
    plot(tout,wheel_torques(:,4),'k-','MarkerSize',6,'LineWidth',2);
    wheel_torques_without_preview = varargin{1};
    plot(tout,wheel_torques_without_preview(:,4),'b--','MarkerSize',6,'LineWidth',2.3);
    legend_text ={'with friction preview','without friction preview'};
    
end

grid on
legend(legend_text,'Location','Best');
% title('Wheels Torque v.s. Time ');
xlabel('Time [s]','fontsize',14)
ylabel('Wheels Torque [Nm]','fontsize',14)
ylim([-1000 1000])
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

