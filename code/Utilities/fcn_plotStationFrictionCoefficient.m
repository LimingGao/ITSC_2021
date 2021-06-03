%%%%%%%%%%%%%%%%%%%%%  Function fcn_plotStationFrictionCoefficient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      plot the tation of cg postion and yaw
%
% Input Variables:
%      P_CG   = the location(ENU) and yaw angle of vehicel center of gravity,format: N*4 matrix
%      friction_coeffs = friction coefficient of four wheels, format: N*4 matrix
%      varargin(1) = estimated friction
%      varargin(2) = tire name 
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

function fcn_plotStationFrictionCoefficient(P_CG,friction_coeffs,varargin)

station_cg = [ 0; cumsum(sqrt(diff(P_CG(:,1)).^2+diff(P_CG(:,2)).^2))]; % calculate station
if ~isempty(varargin) && strcmp(varargin{2},'fl')
    h_fig = figure(847345);
elseif ~isempty(varargin) && strcmp(varargin{2},'rr')
    h_fig = figure(1045);
else
    h_fig = figure(47345);
end
set(h_fig,'Name','fcn_plotStationFrictionCoefficient');
width=600;%
height=400;%
right=100;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
clf;
hold on
if nargin <=3
    plot(station_cg,friction_coeffs(:,1),'ro','MarkerSize',10)
    plot(station_cg,friction_coeffs(:,2),'gs','MarkerSize',6)
    plot(station_cg,friction_coeffs(:,3),'bd','MarkerSize',6)
    plot(station_cg,friction_coeffs(:,4),'cp','MarkerSize',6)
    legend_text = {'Front Left Tire True Friction','Front Right Tire True Friction','Rear Left Tire True Friction','Rear Right Tire True Friction'};
    if nargin ==3
        mu = varargin{1};
        plot(station_cg,mu(:,1),'r.','MarkerSize',10)
        plot(station_cg,mu(:,2),'g.','MarkerSize',10)
        plot(station_cg,mu(:,3),'b.','MarkerSize',10)
        plot(station_cg,mu(:,4),'c.','MarkerSize',10)
        legend_text ={legend_text{1:4}, 'FL estimation','FR estimation','RL estimation','RR estimation'};
        
    end
elseif strcmp(varargin{2},'fl')
    mu = varargin{1};
    plot(station_cg,mu(:,1),'r.','MarkerSize',10)
    plot(station_cg,friction_coeffs(:,1),'b.','MarkerSize',10)
    legend_text = {'Front left tire estimated friction ','Front left tire true friction'};
    
    % Create textarrow
%     annotation(h_fig,'textarrow',[0.286666666666667 0.312],[0.306 0.256],...
%         'String',{'Bridge 1'},'FontName', 'Times New Roman', 'FontSize', 11);
%     annotation(h_fig,'textarrow',[0.453333333333333 0.483333333333335],...
%         [0.868 0.839],'String',{'Bridge 2'},'FontName', 'Times New Roman', 'FontSize', 11);
%     annotation(h_fig,'textarrow',[0.630666666666667 0.659333333333336],...
%         [0.528 0.485],'String',{'Bridge 3'},'FontName', 'Times New Roman', 'FontSize', 11);
elseif strcmp(varargin{2},'rr')
    mu = varargin{1};
    plot(station_cg,mu(:,4),'r.','MarkerSize',10)
    plot(station_cg,friction_coeffs(:,4),'b.','MarkerSize',10)
    legend_text = {'Rear right tire estimated friction ','Rear right tire true friction'};
    
    % Create textarrow
%     annotation(h_fig,'textarrow',[0.286666666666667 0.312],[0.306 0.256],...
%         'String',{'Bridge 1'},'FontName', 'Times New Roman', 'FontSize', 11);
%     annotation(h_fig,'textarrow',[0.453333333333333 0.483333333333335],...
%         [0.868 0.839],'String',{'Bridge 2'},'FontName', 'Times New Roman', 'FontSize', 11);
%     annotation(h_fig,'textarrow',[0.630666666666667 0.659333333333336],...
%         [0.528 0.485],'String',{'Bridge 3'},'FontName', 'Times New Roman', 'FontSize', 11);
else
    error('wrong')
end
legend(legend_text,'Location','best')
% title('Station vs true Friction Coefficient and estimated')
grid on 
xlabel('Station [m]')
ylabel('Friction Coefficient')
xlim([min(station_cg)-10 max(station_cg)+10])
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

