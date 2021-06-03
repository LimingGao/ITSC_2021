%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotTimeCurrentAndPreviewFriction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      plot the time, current friction and preview friction
%
% Input Variables:
%      tout   = time ,format: N*4 matrix
%      friction_coeffs = friction coefficient of four wheels, format: N*4 matrix
%      preview_friction = frictin values in the preveiw area, format: N*5 matrix
%      varargin(1) = tire name
%
% Returned Results:
%
% Example:
% fcn_plotStationCurrentAndPreviewFriction(P_CG,friction_coeffs,preview_friction,varargin)
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

function fcn_plotTimeCurrentAndPreviewFriction(tout,friction_coeffs,preview_friction,varargin)

if ~isempty(varargin) && strcmp(varargin{1},'fl')
    h_fig = figure(3745);
else
    h_fig = figure(3746);
end

set(h_fig,'Name','station of cg postion and friction');
width=600;%
height=420;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on
if nargin == 3
    plot(tout,friction_coeffs(:,1),'ro','MarkerSize',10)
    plot(tout,friction_coeffs(:,2),'gs','MarkerSize',6)
    plot(tout,friction_coeffs(:,3),'yd','MarkerSize',6)
    plot(tout,friction_coeffs(:,4),'cp','MarkerSize',6)
    
    plot(tout,preview_friction(:,3),'b.','MarkerSize',14);
    
    legend('Front Left Tire Current Friction','Front Right Tire Current Friction',...
        'Rear Left Tire Current Friction','Rear Right Tire Current Friction',...
        'Preview Friction','Location','Best','FontSize',10);
    
elseif strcmp(varargin{1},'fl')
    %             plot(tout,preview_friction(:,4),'b.','MarkerSize',10);
    %             plot(tout,preview_friction(:,5),'r.','MarkerSize',8);
    
    plot(tout,friction_coeffs(:,1),'ro','MarkerSize',10)
    %             plot(tout,friction_coeffs(:,2),'gs','MarkerSize',6)
    %             plot(tout,friction_coeffs(:,3),'bd','MarkerSize',6)
    %             plot(tout,friction_coeffs(:,4),'cp','MarkerSize',6)
    
    plot(tout,preview_friction(:,3),'b.','MarkerSize',14);
    
    grid on
    legend('Front Left Tire Current Friction','Preview Friction','Position',[0.478333329280217 0.12598544976951 0.425000008106231 0.0983333354223342], ...
    'FontSize',12)
    
else
    error('wrong')
end

% title('Station vs true Friction Coefficient and estimated')
grid on
xlabel('Time [second]')
ylabel('Friction Coefficient')
% xlim([min(station_cg)-10 max(station_cg)+10])
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

