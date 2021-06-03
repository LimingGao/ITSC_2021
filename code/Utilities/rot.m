%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function rot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   rot returns a 3x3 rotation matrix for rotating a vector about a single 
%   axis. Setting axis = 1 rotates about the e1 axis, axis = 2 rotates 
%   about the e2 axis, axis = 3 rotates about the e3 axis.
% 
% Format:
%   R = rot(angle, axis)
% 
% INPUTS:
%   angle: rotation angle in 'deg'
%   axis: axis of rotation
% 
% OUTPUTS:
%   R: rotation matrix
% 
% Edited by: Satya Prasad, 2019/10/24
% 
% Revision history:
% 2020/06/20
%   - Added input checks. Added figlets to distinguish different
%   sections of the road. Added more comments.
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = rot(angle, axis)
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
    error('rot: Incorrect number of input arguments');
end

if numel(angle) ~= 1
    error('rot: angle must be a scalar in degrees');
end

if numel(axis) ~= 1
    error('rot: axis must be a scalar');
end

%% Compute the rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D0 = pi/180;    % conversion deg to rad
R  = eye(3);    % initialize the rotation matrix
cang = cos(D0*angle);
sang = sin(D0*angle);

% rotation about axis-1 or x-axis
if (1 == axis)
    R(2,2) = cang;
    R(2,3) = sang;
    R(3,2) = -sang;
    R(3,3) = cang;
end
% rotation about axis-2 or y-axis
if (2 == axis)
    R(3,3) = cang;
    R(3,1) = sang;
    R(1,3) = -sang;
    R(1,1) = cang;
end
% rotation about axis-3 or z-axis
if (3 == axis)
    R(1,1) = cang;
    R(1,2) = sang;
    R(2,1) = -sang;
    R(2,2) = cang;
end

end