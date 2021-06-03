%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationFrictionCoeffEstimateNois %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      plot the tation of cg postion and yaw
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      friction_coeffs_estimate_noisy = estimated friction coefficient of four wheels, format: N*4 matrix

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

function fcn_plotStationFrictionCoeffEstimateNoisy(P_CG,friction_coeffs_estimate_noisy)

station_cg = [ 0; cumsum(sqrt(diff(P_CG(:,1)).^2+diff(P_CG(:,2)).^2))]; % calculate station
h_fig = figure(474345);
set(h_fig,'Name','fcn_plotStationFrictionCoeffEstimateNoisy');
width=600;%
height=420;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on
plot(station_cg,friction_coeffs_estimate_noisy(:,1),'r-','MarkerSize',10,'LineWidth',1.5)
plot(station_cg,friction_coeffs_estimate_noisy(:,2),'g--','MarkerSize',5,'LineWidth',1.5)
plot(station_cg,friction_coeffs_estimate_noisy(:,3),'b-','MarkerSize',5,'LineWidth',1.5)
plot(station_cg,friction_coeffs_estimate_noisy(:,4),'c-.','MarkerSize',5,'LineWidth',1.5)
legend_text ={'Front Left Tire','Front Right Tire','Rear Left Tire','Rear Right Tire'};

legend(legend_text,'Location','best')
% title('Station vs true Friction Coefficient and estimated')
grid on
xlabel('Station [m]')
ylabel('Friction Coefficient')
xlim([0 1001])
ylim([0 1])
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
