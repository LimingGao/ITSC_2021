%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotLogitudinalForceWheelSlip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the longitudinal force and wheel slip
%
% Input Variables:
%      Fx   = the longitudinal force of four tires,format: N*4 matrix
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

function fcn_plotLogitudinalForceWheelSlip(Fx,kappa)

h_fig = figure(41115);

set(h_fig,'Name','fcn_plotLogitudinalForceWheelSlip');
width=600;%
height=400;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on

plot(kappa(:,1),Fx(:,1),'g.','MarkerSize',10)
plot(kappa(:,2),Fx(:,2),'r.','MarkerSize',10)
plot(kappa(:,3),Fx(:,3),'co','MarkerSize',10)
plot(kappa(:,4),Fx(:,4),'bo','MarkerSize',10)

legend('front left wheel','front right wheel','rear left wheel','rear right wheel','Location','best')
% title('Station vs true Friction Coefficient and estimated')
grid on
xlabel('wheel longitudinal slip ratio')
ylabel('wheel longitudinal force [N]')

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

