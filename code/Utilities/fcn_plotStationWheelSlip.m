%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationWheelSlip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the station and wheel slip
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      kappa = the wheel slip of four wheels, format: N*4 matrix
% Returned Results:
%
% Example:
% fcn_plotStationFrictionCoefficient(P_CG,friction_coeffs)
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

function fcn_plotStationWheelSlip(P_CG,kappa)

station_cg = [ 0; cumsum(sqrt(diff(P_CG(:,1)).^2+diff(P_CG(:,2)).^2))]; % calculate station

h_fig = figure(41115);

set(h_fig,'Name','fcn_plotTimeWheelSpeed');
width=600;%
height=400;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on

plot(station_cg,kappa(:,1),'g.-','MarkerSize',10)
plot(station_cg,kappa(:,2),'r.','MarkerSize',10)
plot(station_cg,kappa(:,3),'c.','MarkerSize',10)
plot(station_cg,kappa(:,4),'b.-','MarkerSize',10)

legend('front left wheel','front right wheel','rear left wheel','rear right wheel','Location','best')
% title('Station vs true Friction Coefficient and estimated')
grid on
xlabel('station [m]')
ylabel('wheel longitudinal slip')

% axis equal
box on
axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
axes_handle.GridLineStyle = '-.';
axes_handle.GridColor = 'k';
axes_handle.GridAlpha = 0.2;
%             ax.YMinorGrid = 'on';
%             ax.MinorGridAlpha = 0.5
box on
end

