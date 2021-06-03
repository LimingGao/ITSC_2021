function sfcn_enu2lla(block)
% Purpose: coverts cg of the vehicle fron enu to lla
% 
% INPUT 1: 'vehicle_cg' is the cg of the vehicle in ENU
% 
% OUTPUT 1: cg of the vehicle in LLA
% 
% Parameter 1: reference LLA coordinates [-48.876667 -123.393333 0]
% 
% Author: Satya Prasad
% Create Date: 2020-05-21
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
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
% Register the parameters
block.NumDialogPrms     = 1;
block.DialogPrmsTunable = {'Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = 3;

% Override output port properties
block.OutputPort(1).Dimensions = 3;

% Set block sample time to be inherited
block.SampleTimes = [-1, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('Outputs', @Output);   % Required

%% Outputs:
%%   C MEX counterpart: mdlOutputs
function Output(block)
% reference LLAs
ref_lat = block.DialogPrm(1).Data(1);
ref_lon = block.DialogPrm(1).Data(2);
ref_alt = block.DialogPrm(1).Data(3);
% cg of the vehicle in ENU
cg_east  = block.InputPort(1).Data(1);
cg_north = block.InputPort(1).Data(2);
cg_up    = block.InputPort(1).Data(3);

% convert the estimated cg from ENU to LLA
[cg_lat, cg_lon, cg_alt] = enu2geodetic(cg_east, cg_north, cg_up, ref_lat, ref_lon, ref_alt, referenceEllipsoid('wgs84'));

block.OutputPort(1).Data = [cg_lat, cg_lon, cg_alt];
