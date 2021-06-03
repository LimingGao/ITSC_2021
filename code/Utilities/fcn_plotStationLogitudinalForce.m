%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationLogitudinalForce %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the station and vehicle tire longitidinal force 
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      Fx   = the longitudinal force of four tires,format: N*4 matrix
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

function fcn_plotStationLogitudinalForce(P_CG,Fx)

station_cg = [ 0; cumsum(sqrt(diff(P_CG(:,1)).^2+diff(P_CG(:,2)).^2))]; % calculate station
h_fig = figure(41951);

set(h_fig,'Name','fcn_plotStationLogitudinalForce');
width=600;%
height=400;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on

plot(station_cg,Fx(:,1),'g.-','MarkerSize',10)
plot(station_cg,Fx(:,2),'r.','MarkerSize',10)
plot(station_cg,Fx(:,3),'c.','MarkerSize',10)
plot(station_cg,Fx(:,4),'b.-','MarkerSize',10)

legend('front left wheel','front right wheel','rear left wheel','rear right wheel','Location','best')
% title('Station vs true Friction Coefficient and estimated')
grid on
xlabel('Station [m]')
ylabel('Tires longitudianl forces [N] ')

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

