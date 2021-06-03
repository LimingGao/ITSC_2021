%%%%%%%%%%%%%%%%%%% Function fcn_plotSteeringAngle %%%%%%%%%%%%%%%%%%%%%%%%
% fcn_plotSteeringAngle plots steering angle
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
%
% Format:
%   fcn_plotSteeringAngle(delta,tout)
% 
% INPUTS:
%   delta: steering angle in radians
%   tout: simulation time in seconds
% 
% OUTPUTS:
%   plot of steering angle
% 
% Author: Satya Prasad
% Create Date: 2020-06-03
% Revision history:
% 2020/06/03: Added comments and input checks
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcn_plotSteeringAngle(delta,tout)
%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are there right number of inputs?
if nargin ~= 2
    error('Incorrect number of input arguments');
end

% Are both inputs of same size?
if size(delta,1) ~= size(tout,1)
    error('Inputs are of different length');
end
%% plot steering angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_fig = figure(512349);
set(h_fig,'Name','fcn_plotSteeringAngle');
clf;
hold on 
plot(tout,delta(:,1)*180/pi,'b-','LineWidth',1.5);
plot(tout,delta(:,2)*180/pi,'r--','LineWidth',1.5);
plot(tout,delta(:,3)*180/pi,'c-','LineWidth',1.5);
plot(tout,delta(:,4)*180/pi,'k--','LineWidth',1.5);
title('Vehicle Steering Response');
legend('Front Left','Front Right','Rear Left','Rear Right','Location','best');
xlabel('Time (s)');
ylabel('Steering Angle (deg)');

grid on 
%xlim([min(station_cg)-10 1000+10])
axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
axes_handle.GridLineStyle = '-.';
axes_handle.GridColor = 'k';
axes_handle.GridAlpha = 0.2;
%             ax.YMinorGrid = 'on';
%             ax.MinorGridAlpha = 0.5
box on
end