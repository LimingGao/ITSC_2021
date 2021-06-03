%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotLoadFrictionGrid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      vehicle CG trajectory, wheels trajectory, friction grid and friction coefficient at wheels
%
% Input Variables:
%      P_CG,P_FL,P_FR,P_RL,P_RR       position of vehicle CG and four wheels
%      lanesCenter_table
%      frictionGrid_table
%      friction_coeffs        the lane center data , format: table
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

function fcn_plotVehicleTrajectoryAndFrictionGrid(P_CG,P_FL,P_FR,P_RL,P_RR,lanesCenter_table,frictionGrid_table,friction_coeffs)

h_fig = figure(16945);
set(h_fig,'Name','check friction data');
clf;
hold on
%plot3(friction_grid(1:100:end,1),friction_grid(1:100:end,2),friction_grid(1:100:end,4),'g.')
scatter3(frictionGrid_table.position_x(1:11:end),frictionGrid_table.position_y(1:11:end),[frictionGrid_table.friction_coefficient(1:11:end)], [], frictionGrid_table.friction_coefficient(1:11:end),'.');
%plot(lanesCenter_table.center_x,lanesCenter_table.center_y,'k.','MarkerSize',3)
scatter3(lanesCenter_table.center_x,lanesCenter_table.center_y,ones(size(lanesCenter_table.center_x)), 10, ones(size(lanesCenter_table.center_x)),'k.')
%plot(P_CG(:,1),P_CG(:,2),'r','LineWidth',2)
scatter3(P_CG(:,1),P_CG(:,2),ones(size(P_CG(:,1))), 50, ones(size(P_CG(:,1))),'r.')
scatter3(P_FL(:,1),P_FL(:,2),friction_coeffs(:,1), 50, friction_coeffs(:,1),'m.');
scatter3(P_FR(:,1),P_FR(:,2),friction_coeffs(:,2), 50, friction_coeffs(:,2),'m.');
scatter3(P_RL(:,1),P_RL(:,2),friction_coeffs(:,3), 50, friction_coeffs(:,3),'b.');
scatter3(P_RR(:,1),P_RR(:,2),friction_coeffs(:,4), 50, friction_coeffs(:,4),'b.');
% Just for good measure
% colorbar
colormap(winter);
colorbar
legend('friction grid','lane center','vehicle cg trajectory','front left wheel friction','front right wheel friction','rear left wheel friction','rear right wheel friction','Location','best')
grid on
xlabel('xEast')
ylabel('yNorth')
zlabel('friction\_coeficient')
%axis equal
view(0,90)
zlim([0 1])
% axis equal
box on
end

