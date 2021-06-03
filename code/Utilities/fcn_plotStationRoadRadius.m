%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationRoadRadius %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the tation of cg postion and yaw
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      sectionRef_table =  table which has station and curvature columns of
%      road reference curve, format: table
%      road_radius_preview = simulated road radius preview
%      
% Returned Results:
%
% Example:
% fcn_plotStationRoadRadius(P_CG,sectionRef_table,road_radius_preview)
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
function  fcn_plotStationRoadRadius(P_CG,sectionRef_table,road_radius_preview)
h_fig = figure(19982);
set(h_fig,'Name','station and preview road radius');
width=600;%
height=450;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on
station_cg = [0; cumsum(sqrt(sum(diff(P_CG(:,1:2)).^2,2)))]; % calculate station
plot(station_cg,road_radius_preview,'g.','MarkerSize',14);
plot(sectionRef_table.station,1./sectionRef_table.curvature,'r.','MarkerSize',14);
plot(sectionRef_table.station,smoothdata(1./sectionRef_table.curvature,'gaussian',200),'k','MarkerSize',14,'LineWidth',2);

grid on
legend('road\_radius\_preview','road\_radius\_current','road\_radius\_smooth','Position',[0.549222215769026 0.781858043477016 0.355666673119863 0.144572956918397]);
%title('Friction Coefficient in Preview Area v.s. Station ');
xlabel('station [m]','fontsize',14)
ylabel('friction coefficient','fontsize',14)
%ylim([0 1])
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

