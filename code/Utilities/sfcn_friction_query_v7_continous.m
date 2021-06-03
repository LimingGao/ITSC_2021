function sfcn_friction_query_v7(block)
% % Purpose: query the friction values given 4 wheel postion 
% input: friction_grid, road grid and friction, N*4 matrix (E,N,U,friction),meters
%        P_FL, postion of front-left wheel, 1*3 array, (E,N,U), units: meter
%         P_FR,
%         P_RL,
%         P_RR,
% ouput: f, frictions of 4 wheels, 1*4 array[f_FL,f_FR,f_RL,f_RR ]

% Assumption: 2D point search

% Author: Liming Gao
% Create Date: 2020-07-14
% =======update=======
% 1. 
% ======== to do list ============
% 1. 
% 
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
block.NumInputPorts  = 7;
block.NumOutputPorts = 1;

% Register the parameters
block.NumDialogPrms     = 1;
block.DialogPrmsTunable = {'Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = block.DialogPrm(1).Data;
% block.InputPort(1).Dimensions        = 1;
% block.InputPort(1).DatatypeID  = 0;  % double
% block.InputPort(1).Complexity  = 'Real';
% block.InputPort(1).DirectFeedthrough = true;
block.InputPort(2).Dimensions = 4;
block.InputPort(3).Dimensions = 3;
block.InputPort(4).Dimensions = 3;
block.InputPort(5).Dimensions = 3;
block.InputPort(6).Dimensions = 3;
block.InputPort(7).Dimensions = 1;

% Override output port properties
block.OutputPort(1).Dimensions = 4;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
% Set block sample time to be inherited
block.SampleTimes = [0, 0];

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
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);  
block.RegBlockMethod('InitializeConditions', @InitConditions);  
block.RegBlockMethod('Outputs', @Output);   % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C MEX counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  %block.Dwork(1).UsedAsDiscState = true;
  
%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C MEX counterpart: mdlInitializeConditions
%%
function InitConditions(block)

  %% Initialize Dwork
  block.Dwork(1).Data = 1.0;

%% Outputs: this function is executed in each simluation step
%%   C MEX counterpart: mdlOutputs
function Output(block)
% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
% variable

persistent flag; 
persistent lane_grid_cluster;
persistent friction_cluster;
persistent friction_cluster_range;
persistent vehi_wheelP_query;
persistent friction_wheels;
persistent Md_lane_grid;  % MD of all data
persistent Md_NS; % MD of range data
persistent Idx;
persistent Distance; % find all the grid in this distance around the vehicle
persistent range_overlap; % the distance vechicle move for next surrounding range query
persistent counts; % station from the start point 
persistent P_CG_last; % slast position of the vehicle
persistent Road_section_last; % slast position of the vehicle
persistent CG_station; % station from the start point 

% set the flag to empty if there is a new section friction
% if ~isempty(lane_grid_cluster)
if Road_section_last ~= block.InputPort(7).Data
    flag = [];
end
% end
Road_section_last = block.InputPort(7).Data; % 

vehi_wheelP_query =[block.InputPort(3).Data(1:2) block.InputPort(4).Data(1:2) block.InputPort(5).Data(1:2) block.InputPort(6).Data(1:2)]'; % 4*2 matrix

% Initialization
if isempty(flag)
    flag = 1;
    
    counts =1;
    Distance = 100; %meters, radius
    range_overlap = 60;  % meters
    
    % step 1: prepare data
    lane_grid_cluster = block.InputPort(1).Data(:,1:2); % xEast and yNorth, ignore zUp
    friction_cluster = block.InputPort(1).Data(:,4);  % friction coefficient value
    Md_lane_grid= KDTreeSearcher(lane_grid_cluster); % prepare data tree for range search
    
    % step2: first search. Query the surrounding range
    Idx_range = rangesearch(Md_lane_grid,block.InputPort(2).Data(1:2)',Distance,'SortIndices',false); %For faster computation, specify to keep the indices of the nearest neighbors unsorted.
    
    Md_NS = createns(lane_grid_cluster(Idx_range{1},:)); % Create nearest neighbor searcher object
    friction_cluster_range = friction_cluster(Idx_range{1});
    
    %nb_grid = block.InputPort(1).Dimensions(1);
    CG_station = 0;
    P_CG_last = block.InputPort(2).Data(1:2)';
end

CG_station = CG_station + sqrt(sum((block.InputPort(2).Data(1:2)' - P_CG_last).^2));  % % just use E and N

P_CG_last = block.InputPort(2).Data(1:2)'; 

if CG_station > counts*range_overlap  % update search range
    counts = counts + 1;
    Idx_range = rangesearch(Md_lane_grid,block.InputPort(2).Data(1:2)',Distance,'SortIndices',false); %For faster computation, specify to keep the indices of the nearest neighbors unsorted.
    Md_NS = createns(lane_grid_cluster(Idx_range{1},:)); % Create nearest neighbor searcher object
    friction_cluster_range = friction_cluster(Idx_range{1});
end

% step3: query wheel nearest grid just in the surrounding range
Idx = knnsearch(Md_NS,vehi_wheelP_query);

% step4: find the nearest friction values for each wheel position
friction_wheels = friction_cluster_range(Idx);

if isempty(friction_wheels)
    friction_wheels = zeros(4,1);
    %warning('Wheel Friction can not be found!')
end

block.OutputPort(1).Data = friction_wheels;