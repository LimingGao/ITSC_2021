%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotLoadFrictionGrid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      check queried/loaded friction grid and lane center line
%
% Input Variables:
%      frictionGrid_table       = the friction grid data, format: table
%      lanesCenter_table         = the lane center data , format: table
%      frictionGrid_Cluster_table   = the clusterd friction grid data, format: table
% Returned Results:
%
% Example:
% fcn_plotLoadFrictionGrid(frictionGrid_table,results_lane_node)
% 
% Processing Flow:
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

function fcn_plotLoadFrictionGrid(frictionGrid_table,lanesCenter_table,frictionGrid_Cluster_table)

h_fig = figure(458745);
set(h_fig,'Name','check friction grid and lane center line');
clf;
hold on
%plot(frictionGrid_table.position_x,frictionGrid_table.position_y,'g.','MarkerSize',15)
scatter3(frictionGrid_table.position_x,frictionGrid_table.position_y,frictionGrid_table.friction_coefficient, [], frictionGrid_table.friction_coefficient,'.');
%plot(lanesCenter_table.center_x,lanesCenter_table.center_y,'k.','MarkerSize',3)
scatter3(lanesCenter_table.center_x,lanesCenter_table.center_y,ones(size(lanesCenter_table.center_x)),18, ones(size(lanesCenter_table.center_x)),'k.')

if nargin > 2
    %plot(frictionGrid_Cluster_table.position_x,frictionGrid_Cluster_table.position_y,'r.','MarkerSize',3)
    scatter3(frictionGrid_Cluster_table.position_x,frictionGrid_Cluster_table.position_y,ones(size(frictionGrid_Cluster_table.position_x)),18,'r.')
    legend('friction grid','lane center','friction grid cluster')
else
    legend('friction grid','lane center')
end
colormap default
grid on
xlabel('xEast')
ylabel('yNorth')
zlabel('friction\_coeficient -1 ')
%axis equal
view(0,90)
zlim([0 1])
% axis equal
box on
hold off
end

