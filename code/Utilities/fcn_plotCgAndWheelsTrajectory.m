%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotLoadFrictionGrid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      check queried/loaded friction grid and lane center line
%
% Input Variables:
%      frictionGrid_table       the friction grid data, format: table
%      results_lane_node        the lane center data , format: table
% 
% Returned Results:
%
% Example:
% fcn_plotLoadFrictionGrid(frictionGrid_table,results_lane_node)
% 
% Processing Flow:
%
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

function fcn_plotCgAndWheelsTrajectory(frictionGrid_table,lanesCenter_table,P_CG,P_FL,P_FR,P_RL,P_RR)

h_fig = figure(4745);
set(h_fig,'Name','check data trajectory and wheel position');
clf;
hold on
plot(frictionGrid_table.position_x(1:11:end),frictionGrid_table.position_y(1:11:end),'g.','MarkerSize',15)
plot(lanesCenter_table.center_x,lanesCenter_table.center_y,'k.','MarkerSize',3)

plot(P_CG(:,1),P_CG(:,2),'r','LineWidth',1)

plot(P_FL(:,1),P_FL(:,2),'m.','MarkerSize',10)
plot(P_FR(:,1),P_FR(:,2),'m.','MarkerSize',10)
plot(P_RL(:,1),P_RL(:,2),'b.','MarkerSize',5)
plot(P_RR(:,1),P_RR(:,2),'b.','MarkerSize',5)
legend('friction grid','lane center','vehicle cg trajectory','front left','front right','rear left','rear right')
grid on
xlabel('xEast [m]')
ylabel('yNorth [m]')
% xlim([min(P_CG(:,1))-10 max(P_CG(:,1))+10])
% axis equal
box on
end

