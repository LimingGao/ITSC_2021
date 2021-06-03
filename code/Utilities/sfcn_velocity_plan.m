% % Purpose: plan velocity profile from current postion to preview position
%            with constant acceleration assumption
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

function sfcn_velocity_plan(block)

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
persistent Vc_last; % last command speed
persistent preview_time; % time duration of the planner
persistent horizon; % how many steps the velocity plan
persistent V_horizon; % Velocity at the horizon

Vx = block.InputPort(1).Data;
Va_ahead = block.InputPort(2).Data;
Vd = block.InputPort(3).Data;

% Initialization
if isempty(flag)
    flag = 1;
    
    preview_distance = 100; %meters, radius
    sample_time = 0.01; % sample time of 100Hz
    % step 1: 
    if Vd > Va_ahead
        % calculate the constant acceleration from current postion to preview position
        acceleration = (Va_ahead-Vx)*(Va_ahead+Vx)/(2*preview_distance);
        
        % calcualte the horizon
        preview_time = round((2*preview_distance)/(Va_ahead+Vx),2);
        horizon = preview_time/sample_time;
        V_horizon = 111;
        
        % calculate the speed of next time step
        Vc = Vx + acceleration*sample_time;
    else
        Vc = Vd;
    end
    
    Vc_last = Vc;
else
    % step 1:
    if Vd > Va_ahead
        % calculate the constant acceleration from current postion to preview position
        acceleration = (Va_ahead-Vx)*(Va_ahead+Vx)/(2*preview_distance);
        % calculate the speed of next time step
        Vc = Vx + acceleration*sample_time;
    else
        Vc = Vd;
    end
    % step 2: choose the smaller one
    if Vc > Vc_last
        Vc = Vc_last;
    end
    Vc_last = Vc;

end

% if isempty(friction_preview)
%     friction_preview = 0;
%     warning('Wheel Friction can not be found!')
% end

block.OutputPort(1).Data = Vc;