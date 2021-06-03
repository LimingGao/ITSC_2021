%%%%%%%%%%%%%%%%%%% Function fcn_plotLateralOffsetError %%%%%%%%%%%%%%%%%%%
% fcn_plotLateralOffsetError plots lateral offset error
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
%
% Format:
%   fcn_plotLateralOffsetError(errorLateralOffset,tout)
% 
% INPUTS:
%   errorLateralOffset: lateral offset error
%   tout: simulation time
% 
% OUTPUTS:
%   plot of lateral offset error
% 
% Author: Satya Prasad
% Create Date: 2020-06-03
% Revision history:
% 2020/06/03: Added comments and input checks
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcn_plotspeed_cmd_actual(speed,speed_without_preview,station,errorLateralOffset_without_preview);
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

% Are both inputs of same size?
if length(speed_without_preview) ~= length(station)
    error('Inputs are of different length');
end
%% plot lateral offset error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fig = figure(5239768);
width=600;%
height=250;%
right=150;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
set(h_fig,'Name','fcn_plotLateralOffsetError_both');
clf;
hold on 
index = find(abs(errorLateralOffset_without_preview) > 2.5);
offset = 500;
plot(station,speed,'k-','MarkerSize',6,'LineWidth',1.5);
plot(station,speed_without_preview,'b:','MarkerSize',6,'LineWidth',2.3);


annotation(h_fig,'textbox',...
    [0.790666666666669 0.236179487179487 0.082666666666666 0.120192307692308],...
    'String',{'(b)'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
% title('Lateral Offset Error');
legend('actual speed ','command speed',...
    'Position',[0.547416663408281 0.748404592559685 0.358333339850108 0.178285260231067])
xlabel('Station (m)');
ylabel('Speed (m/s)');
ylim([0 30])
xlim([0 2100])
grid on
axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
axes_handle.GridLineStyle = '-.';
axes_handle.GridColor = 'k';
axes_handle.GridAlpha = 0.2;
%             ax.YMinorGrid = 'on';
%             ax.MinorGridAlpha = 0.5
box on
end