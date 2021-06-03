%%%%%%%%%%%%%%%%%%% Function fcn_plotLateralOffsetError %%%%%%%%%%%%%%%%%%%
% fcn_plotLateralOffsetError plots lateral offset error
%
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
%
% Format:
%   fcn_plotLateralOffsetError(errorLateralOffset,tout)
%
% INPUTS:
%   errorLateralOffset: look ahead lateral offset error
%   tout: simulation time
%   lateral_error: current lateral error
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
function fcn_plotLateralOffsetError(tout,errorLateralOffset,varargin)
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

% Are both inputs of same size?
if length(errorLateralOffset) ~= length(tout)
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
h_fig = figure(5123648);
width=600;%
height=450;%
right=150;%
bottom=400;%
set(gcf,'position',[right,bottom,width,height])
set(h_fig,'Name','Lateral Offset Error');
clf;
hold on
if nargin ==3
    lateral_error = varargin{1};
    plot(tout,errorLateralOffset,'b-','MarkerSize',6,'LineWidth',1.5);
    plot(tout,lateral_error,'k-','MarkerSize',6,'LineWidth',1.5);
    %title('Lateral Lookahead Offset Error');
    legend('Lateral Look Ahead Error','Actual Lateral Error','Location','best');
elseif nargin ==2
    plot(tout,errorLateralOffset,'b-','MarkerSize',6,'LineWidth',1.5);
    %title('Lateral Lookahead Offset Error');
    legend('Lateral Look Ahead Error','Location','best');
else
    error('Incorrect number of input arguments');
end
xlabel('Time (s)');
ylabel('Error (m)');
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