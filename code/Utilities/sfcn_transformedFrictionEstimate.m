%%%%%%%%%%%%%%% S-Function sfcn_transformedFrictionEstimate %%%%%%%%%%%%%%%%
% Purpose:
%   Estimates the transformed friction
% 
% INPUTS:
%   INPUT 1: friction value from the DB sensor
%   Input 2: intercept and slope of the transformation estimate
% 
% OUTPUTS:
%   OUTPUT 1: transformed friction estimate
% 
% Author: Satya Prasad
% Create Date: 2020-10-01
% Revision history:
% 2020\10\01: 
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_transformedFrictionEstimate(block)

% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 1;

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = 4;
block.InputPort(2).Dimensions = [4, 2];

% Override output port properties
block.OutputPort(1).Dimensions = 4;

% Set block sample time to be inherited
block.SampleTimes = [-1, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('Outputs', @Output); % Required

function Output(block)
%% Output function to estimate transformed friction

mu_FL = block.InputPort(2).Data(1,2)*block.InputPort(1).Data(1) + block.InputPort(2).Data(1,1);
mu_FR = block.InputPort(2).Data(2,2)*block.InputPort(1).Data(2) + block.InputPort(2).Data(2,1);
mu_RL = block.InputPort(2).Data(3,2)*block.InputPort(1).Data(3) + block.InputPort(2).Data(3,1);
mu_RR = block.InputPort(2).Data(4,2)*block.InputPort(1).Data(4) + block.InputPort(2).Data(4,1);

block.OutputPort(1).Data = [mu_FL; mu_FR; mu_RL; mu_RR];