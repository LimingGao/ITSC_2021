% % Purpose: plan velocity profile from current postion to preview position
%            with constant acceleration assumption(space based)
% 
% input: Vx: current vehicle velocity, scalar,m/s
%        Va_ahead: allowable maximum velocity ahead of vehicle
%        Vd: desired velocty to track the trajectory 
%        
% ouput: Vc, command velocity to the P controller

% Assumption: 2D point search

% Author: Liming Gao
% Create Date: 2020-09-10
% =======update=======
% 1. 
% ======== to do list ============
% 1.
% 

function sfcn_velocity_plan_v2(block)

% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C MEX counterpart: mdlInitializeSizes
%%   It seems that this function is executed once at the start of simluation
function setup(block) 
% Register number of input and output ports
block.NumInputPorts  = 3;
block.NumOutputPorts = 1;

% Register the parameters
% block.NumDialogPrms     = 2;
% block.DialogPrmsTunable = {'Nontunable','Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = 1;
% block.InputPort(1).Dimensions        = 1;
% block.InputPort(1).DatatypeID  = 0;  % double
% block.InputPort(1).Complexity  = 'Real';
% block.InputPort(1).DirectFeedthrough = true;
block.InputPort(2).Dimensions = 1;
block.InputPort(3).Dimensions = 1;

% Override output port properties
block.OutputPort(1).Dimensions = 1;
block.OutputPort(1).SamplingMode = 'Sample';
% block.OutputPort(2).Dimensions = 2;
% block.OutputPort(2).SamplingMode = 'Sample';
% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
% Set block sample time to be inherited
block.SampleTimes = [-1, 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

%%-----------------------------------------------------------------
%%The MATLAB S-function uses an internal registry for all
%%block methods. You should register all relevant methods
%%(optional and required) as illustrated below. You may choose
%%any suitable name for the methods and implement these methods
%%as local functions within the same file. See comments
%%provided for each function for more information.
%%-----------------------------------------------------------------
% Register methods
block.RegBlockMethod('Outputs', @Output);   % Required
%%%%%%%%%%%%%%%%%%%% need this method if theer are more than one outpuuts
% SetInputPortSamplingMode:
%   Functionality    : Check and set input and output port
%                      attributes and specify whether the port is operating
%                      in sample-based or frame-based mode
%   C MEX counterpart: mdlSetInputPortFrameData.
%   (The DSP System Toolbox is required to set a port as frame-based)
%
% block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);

%end setup
%% Outputs: this function is executed in each simluation step
%%   C MEX counterpart: mdlOutputs
function Output(block)
% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
% variable

persistent flag; 
% input variable 
persistent Vx; % current vehicle velocity, scalar,m/s
persistent Va_ahead;% allowable maximum velocity ahead of vehicle
persistent Vd; %desired velocty to track the trajectory

% output variable 
persistent Vc; % command velocity to the P controller

% intermediate variable
persistent sample_time; % the smaple time of the simulink 
persistent acceleration; % the acceleration from current postion to preview position
persistent preview_distance; % the distance from vehicle cg point to preview point

% persistent preview_time; % time duration of the planner
% persistent horizon; % how many steps the velocity plan
% persistent V_horizon; % Velocity at the horizon
persistent station; % station 
persistent station_delta; % station 
persistent flag_brake;

persistent Va_ahead_buffer; 
persistent Vx_buffer;
persistent buffer_index; 
persistent counts; 

Vx = block.InputPort(1).Data;
Va_ahead = block.InputPort(2).Data;
Vd = block.InputPort(3).Data;

% Initialization
if isempty(flag)
    flag = 1;
    counts = 1;
    preview_distance = 120; %meters, radius
    sample_time = 0.01; % sample time of 100Hz
    station = 0; % initlal

    flag_brake = 0;
    buffer_index = 1;
    station_delta = 0;
    Va_ahead_buffer= [];
    Vx_buffer= [];
end

counts = counts +1;
% if mod(counts,100)==0
%     warning('I am ok!')
% end

% step 1:
if Vd > Va_ahead
    flag_brake =1;
    Vx_buffer = [Vx_buffer Vx];
    Va_ahead_buffer = [Va_ahead_buffer Va_ahead];
    station_delta = station_delta + Vx*sample_time;
else
    
end

if flag_brake ==1
    station = station + Vx*sample_time;
end

if station > 0 && station <= preview_distance
    
    % calculate the constant acceleration from current postion to preview position
    acceleration = (min(Va_ahead_buffer)-max(Vx_buffer))*(min(Va_ahead_buffer)+max(Vx_buffer))/(2*preview_distance);
    % calculate the speed of next time step
    Vc = Vx + 35*acceleration*sample_time;
elseif station > preview_distance && station <= (preview_distance + 15+ station_delta)
    %Vc = Va_ahead_buffer(buffer_index);
    
    Vc = min(Va_ahead_buffer); % more convservation 
%     if counts > 2000
%         Vc = min(Va_ahead_buffer)-1; % more convservation 
%     end
    buffer_index = min(buffer_index + 1,length(Va_ahead_buffer));
    
elseif station > (preview_distance + 15+ station_delta)
    station_delta = 0;
    Va_ahead_buffer = [];
    Vx_buffer= [];
    flag_brake = 0;
    station = 0 ;
    buffer_index =1;
end

if flag_brake ==0
    
    Vc = Vd;
end

% if isempty(friction_preview)
%     friction_preview = 0;
%     warning('Wheel Friction can not be found!')
% end

block.OutputPort(1).Data = Vc;