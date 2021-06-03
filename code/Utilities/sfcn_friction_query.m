function sfcn_friction_query(block)
% % Purpose: query the friction values given 4 wheel postion 
% input: friction_grid, road grid and friction, N*4 matrix (E,N,U,friction),meters
%        P_FL, postion of front-left wheel, 1*3 array, (E,N,U), units: meter
%         P_FR,
%         P_RL,
%         P_RR,
% ouput: f, frictions of 4 wheels, 1*4 array[f_FL,f_FR,f_RL,f_RR ]

% Assumption: 2D point search

% Author: Liming Gao
% Create Date: 2020-05-20
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

function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 5;
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
block.InputPort(2).Dimensions = 3;
block.InputPort(3).Dimensions = 3;
block.InputPort(4).Dimensions = 3;
block.InputPort(5).Dimensions = 3;

% Override output port properties
block.OutputPort(1).Dimensions = 4;

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

%end setup
%% Outputs:
%%   C MEX counterpart: mdlOutputs
function Output(block)
% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
% variable

persistent flag; 
persistent lane_grid_sorted;
persistent friction_sorted;
persistent vehi_wheelP_query;
persistent friction_wheels;
persistent Idx;
persistent Distance; % nearest distance
persistent next_search_index; % search index range
persistent start_index; % first element of search index range
persistent nb_grid;

vehi_wheelP_query =[block.InputPort(2).Data(1:2) block.InputPort(3).Data(1:2) block.InputPort(4).Data(1:2) block.InputPort(5).Data(1:2)]';

% Initialization
if isempty(flag)
    flag = 1;
    Idx = zeros(4,1);
    Distance = zeros(4,1);
    next_search_index = zeros(1000,1);
    
    %nb_grid = length(block.InputPort(1).Data);
    nb_grid = block.InputPort(1).Dimensions(1);
    lane_grid_sorted = block.InputPort(1).Data(:,1:2); % xEast and yNorth, ignore zUp
    friction_sorted = block.InputPort(1).Data(:,4); 
 
    [Idx,Distance] = knnsearch(lane_grid_sorted,vehi_wheelP_query); % resolve the error: mxArray indices are not supported. Consider first storing the index expression into a non-mxarray temporary.
    
    friction_wheels = block.InputPort(1).Data(Idx,4);
    next_search_index = max(1,min(Idx)-500):min(max(Idx) + 8000,nb_grid);
    
else
    
    % step1: find the nearest grid for each wheel

    [Idx,Distance] = knnsearch(lane_grid_sorted(next_search_index,:),vehi_wheelP_query);
    
    % Notes: query error check, but it needs furtehr 
    %     if any(D > 0.5)
    %         %error('nearest point can not be found')
    %         fprintf(1,'step\n');
    %     end
    
    
    % step2: find the nearest friction values for each wheel position
    friction_wheels = friction_sorted(next_search_index(Idx));

    % step3: update the next searching range
    start_index = next_search_index(1);
    next_search_index = max(1,min(Idx)+ start_index - 500):min(max(Idx)+start_index + 8000,nb_grid);

end

block.OutputPort(1).Data = friction_wheels;

