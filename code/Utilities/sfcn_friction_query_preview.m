% % Purpose: query the ahead friction values
% input: friction_grid: road grid and friction, N*4 matrix (E,N,U,friction),meters
%        ref_trajectory: pre planned reference trajectory, N*20 matrix (E,N,U,station,...),meters
%        P_CG: postion of vehicle center, 1*3 array, (E,N,U), units: meter
%        heading: yaw angle of vehicle
%        :
%        :
% ouput: friction_preview, scalar

% Assumption: 2D point search

% Author: Liming Gao
% Create Date: 2020-09-07
% =======update=======
% 1. 
% ======== to do list ============
% 1. Make the preview_distance to be adaptive
% 2. for the nearest station query,using pojection for higher accuracy
% 3. how to represent the preview friction


function sfcn_friction_query_preview(block)

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
block.NumInputPorts  = 5;
block.NumOutputPorts = 1;

% Register the parameters
block.NumDialogPrms     = 3;
block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = block.DialogPrm(1).Data;
% block.InputPort(1).Dimensions        = 1;
% block.InputPort(1).DatatypeID  = 0;  % double
% block.InputPort(1).Complexity  = 'Real';
% block.InputPort(1).DirectFeedthrough = true;
block.InputPort(2).Dimensions = block.DialogPrm(2).Data;
block.InputPort(3).Dimensions = 4;
block.InputPort(4).Dimensions = 1;
block.InputPort(5).Dimensions = 1;

% Override output port properties
block.OutputPort(1).Dimensions = 5;
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
persistent lane_grid_cluster; % the coodinate of grid (xEast,yNorth)
persistent friction_cluster; % the friction values at each grid (friction_coefficient)
persistent friction_cluster_range;
persistent Md_lane_grid;  % MD of all data
persistent Md_NS; % MD of range data

persistent ref_trajectory; % the coodinate of reference trajectory(xEast,yNorth)
persistent ref_trajectory_station; % the station of reference trajectory (station)
persistent Md_ref_trajectory; % MD of reference trajectory
persistent Idx_station; % index of station

persistent Current_station; % vehicle current station at the reference trajectory
persistent preview_distance; % the preview distance for friction query
persistent Preview_station; % station at the preview location
persistent P_preview; % the coodinate of the preview point
persistent preview_range; % the radius of preview area 
persistent friction_preview; % the value to represent the preview friction 
persistent friction_preview_area; % the friction values in the preview area 

persistent range_distance; % find all the data in this distance around the preview point
persistent range_overlap; % the distance vechicle move for next surrounding range query
persistent counts; % station from the start point 
persistent P_CG_last; % last position of the vehicle
persistent P_CG; % current position of the vehicle
persistent Road_section_last; % slast position of the vehicle
persistent CG_station; % station from the start point 

persistent Idx_preview_grid; % the index of gird to represent the preview area

% set the flag to empty if there is a new section friction
% if ~isempty(lane_grid_cluster)
if Road_section_last ~= block.InputPort(4).Data
    flag = [];
end
% end
Road_section_last = block.InputPort(4).Data; % 

Desired_station = block.InputPort(5).Data;

% Initialization
if isempty(flag)
    flag = 1;
    
    counts =1;
    range_distance = 100; %meters, radius
    range_overlap = 50;  % meters
    preview_distance = block.DialogPrm(3).Data; %meters
    preview_range = 1.5; % the radius of preview circle range, meters
    
    % step 1: prepare friction map data
    lane_grid_cluster = block.InputPort(1).Data(:,1:2); % xEast and yNorth, ignore zUp
    friction_cluster = block.InputPort(1).Data(:,4);  % friction coefficient value
    Md_lane_grid= KDTreeSearcher(lane_grid_cluster); % prepare data tree for range search
    
    % step 2: prepare reference trajectory data
    ref_trajectory = block.InputPort(2).Data(:,12:13); % xEast and yNorth, ignore zUp
    ref_trajectory_station = block.InputPort(2).Data(:,18);  % total station of the refrence trajectory, it is a preplaning trajectory
    Md_ref_trajectory= KDTreeSearcher(ref_trajectory); % prepare data tree for trajectory search 
     
    % step3: first search. query the nearest station values
    Idx_station = knnsearch(Md_ref_trajectory,block.InputPort(3).Data(1:2)'); %%current vehicle postion: block.InputPort(3).Data(1:2)';
    Current_station = ref_trajectory_station(Idx_station); % Note: this should be done using projection for higher accuracy
    
    % step4: find the station and coordinate at the preview point
    Preview_station = Current_station + preview_distance; 
    [~,Idx_preview_station] = mink(abs(Preview_station-ref_trajectory_station),1);
    P_preview = ref_trajectory(Idx_preview_station,:);
    
    % step5: determine the search range centered at the preview point
    Idx_range = rangesearch(Md_lane_grid,P_preview,range_distance,'SortIndices',false); %For faster computation, specify to keep the indices of the nearest neighbors unsorted.
    
    Md_NS = createns(lane_grid_cluster(Idx_range{1},:)); % Create nearest neighbor searcher object
    friction_cluster_range = friction_cluster(Idx_range{1});
    
    CG_station = 0;
    P_CG_last = block.InputPort(3).Data(1:2)';
end

% step3:query the nearest station values to find current station 
Idx_station = knnsearch(Md_ref_trajectory,block.InputPort(3).Data(1:2)');
Current_station = ref_trajectory_station(Idx_station); % Note: this should be done using projection for higher accuracy

% step4: find the station and coordinate at the preview point
Preview_station = Current_station + preview_distance;
[~,Idx_preview_station] = mink(abs(Preview_station-ref_trajectory_station),1);
P_preview = ref_trajectory(Idx_preview_station,:);

% vehicle station: the distance vehicle has traveled
CG_station = CG_station + sqrt(sum((block.InputPort(3).Data(1:2)' - P_CG_last).^2));  % % just use E and N

P_CG_last = block.InputPort(3).Data(1:2)'; 

% update friction search range once the distance of vehicle traveled pass a threshold, wheraas we use global serach for trajectory
if CG_station > counts*range_overlap  
    counts = counts + 1;
    %repeat steps 5 to update search range 
    % step5: range search centered at the preview point
    Idx_range = rangesearch(Md_lane_grid,P_preview,range_distance,'SortIndices',false); %For faster computation, specify to keep the indices of the nearest neighbors unsorted.
    Md_NS = createns(lane_grid_cluster(Idx_range{1},:)); % Create nearest neighbor searcher object
    friction_cluster_range = friction_cluster(Idx_range{1});
end

% step6: query preview grids area in surrounding range
% Idx = knnsearch(Md_NS,vehi_wheelP_query);
Idx_preview_grid = rangesearch(Md_NS,P_preview,preview_range,'SortIndices',false); %For faster computation, specify to keep the indices of the nearest neighbors unsorted.
    
% step7: find the friction in the preview area
friction_preview_area = friction_cluster_range(Idx_preview_grid{1});
friction_preview = [min(friction_preview_area) mean(friction_preview_area) max(friction_preview_area)];

if isempty(friction_preview) || any(isnan(friction_preview))
    friction_preview = [0 0 0];
    %warning('Wheel Friction can not be found!')
end

block.OutputPort(1).Data = [P_preview friction_preview];