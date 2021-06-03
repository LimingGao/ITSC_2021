%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationFrictionPreview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      plot the tation of cg postion and preview friction
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      preview_friction = preview friction ahead of vehicle, format: N*5 matrix
%      
% Returned Results:
%
% Example:
%  fcn_plotStationFrictionPreview(P_CG,preview_friction)
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

function fcn_plotStationFrictionPreview(P_CG,preview_friction)

h_fig = figure(3782);
set(h_fig,'Name','station and preview friction');
width=600;%
height=450;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on
station_cg = [0; cumsum(sqrt(sum(diff(P_CG(:,1:2)).^2,2)))]; % calculate station
plot(station_cg,preview_friction(:,3),'g.','MarkerSize',14);
plot(station_cg,preview_friction(:,4),'b.','MarkerSize',10);
plot(station_cg,preview_friction(:,5),'r.','MarkerSize',8);
grid on
legend('Minmum Friction Coeff','Mean Friction Coeff', 'Max Friction Coeff','Position',[0.549222215769026 0.781858043477016 0.355666673119863 0.144572956918397]);
title('Friction Coefficient in Preview Area v.s. Station ');
xlabel('station [m]','fontsize',14)
ylabel('friction coefficient','fontsize',14)
ylim([0 1])
xlim([0 1050])
axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
axes_handle.GridLineStyle = '-.';
axes_handle.GridColor = 'k';
axes_handle.GridAlpha = 0.2;
%             ax.YMinorGrid = 'on';
%             ax.MinorGridAlpha = 0.5
box on

end

