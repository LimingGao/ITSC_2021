%%%%%%%%%%%%%% Function fcn_checkValidityUnixAndAimsunRange %%%%%%%%%%%%%%%
% Purpose:
%   fcn_checkValidityUnixAndAimsunRange checks the validity of unix time,
%   lower and upper bounds on system entrance time by comparing structure
%   'time' with 'unixTimeAndAimsunTimeRange'
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   fcn_checkValidityUnixAndAimsunRange(unixTimeAndAimsunTimeRange,time)
% 
% INPUTS:
%   unixTimeAndAimsunTimeRange: It's a Nx3 vector. First column contains
%   Unix Time for trips on 'names.road'. Second and third columns contains 
%   lower and upper bounds on System Entrance Time respectively for trips 
%   on 'names.road'
%   time: It's a structure containing unix time, lower and upper bounds on
%   system entrance time
% 
% 
% Author: Satya Prasad
% Create Date: 2020-06-03
% Revision history:
% 
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcn_checkValidityUnixAndAimsunRange(unixTimeAndAimsunTimeRange,time)
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

%% Check the validity of Unix Time and Bounds on System Entrance Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the existence of given unix time in the valid set of unix time(s)
if ~ismember(time.Unix,unixTimeAndAimsunTimeRange(:,1))
    error('ERROR: Given UNIX time is not in the valid set of UNIX time(s)');
else
    fprintf(1,'\nGiven UNIX time = %d is in the valid set of UNIX time(s) \n', time.Unix);
    
    index = find(unixTimeAndAimsunTimeRange(:,1),time.Unix);
    % Check the validity of given global time range
    if unixTimeAndAimsunTimeRange(index,2) <= time.AimsunLB || ...
            unixTimeAndAimsunTimeRange(index,3) >= time.AimsunLB
        fprintf(1,'time.AimsunLB = %.2f is within the bounds on system entrance time \n',time.AimsunLB);
        
        if unixTimeAndAimsunTimeRange(index,2) <= time.AimsunUB || ...
            unixTimeAndAimsunTimeRange(index,3) >= time.AimsunUB
            fprintf(1,'time.AimsunUB = %.2f is within the bounds on system entrance time \n',time.AimsunUB);
        else
            error('ERROR: time.AimsunLB is not within the bounds on system entrance time');
        end
        
    else
        error('ERROR: time.AimsunLB is not within the bounds on system entrance time');
    end
end

end