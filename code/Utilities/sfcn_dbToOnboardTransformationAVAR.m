%%%%%%%%%%%%%% S-Function sfcn_dbToOnboardTransformationAVAR %%%%%%%%%%%%%%%
% Purpose:
%   Estimates the transformation to transform friction value from DB sensor
%   to a value from Onboard sensor
% 
% INPUTS:
%   INPUT 1: friction value from the DB sensor
%   Input 2: friction value from the Onboard sensor
%   Input 3: window size for estimating the transformation
% 
% OUTPUTS:
%   OUTPUT 1: intercept and slope of the transformation estimate
% 
% PARAMETERS:
%   PARAMETER 1: horizon length - decides the maximum length of the data to
%   be stored.
% 
% Author: Satya Prasad
% Create Date: 2020-10-06
% Revision history:
% 2020\10\02: 
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_dbToOnboardTransformationAVAR(block)

% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

function setup(block)
%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
% Register number of input and output ports
block.NumInputPorts  = 3;
block.NumOutputPorts = 1;

% Register the parameters
block.NumDialogPrms     = 1;
block.DialogPrmsTunable = {'Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = 4;
block.InputPort(2).Dimensions = 4;
block.InputPort(3).Dimensions = 1;

% Override output port properties
block.OutputPort(1).Dimensions = [4, 2];

% Set block sample time to be inherited
block.SampleTimes = [-1, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitConditions);
block.RegBlockMethod('Outputs', @Output); % Required

function DoPostPropSetup(block)
%% Setup Dwork
block.NumDworks = 1;
block.Dwork(1).Name = 'flag'; % variable to do the initialization
block.Dwork(1).Dimensions = 1;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = 'Real';

function InitConditions(block)
%% Initialize Dwork
block.Dwork(1).Data = 1; % acts as a flag

function Output(block)
%% Output function to estimate transformation to transform friction from DB TO Onboard

% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
% variables to store friction values measured by the Onboard sensor
persistent mu_onboard;
% variables to store friction values from the DB sensor
persistent mu_db
% % variables to store intercept and slope of transformation
% persistent transformation

if block.Dwork(1).Data
    % Intialize the variables for every new simulation run
    mu_db = [];
    mu_onboard = [];
    % reset the variable
    block.Dwork(1).Data = 0;
end

% window size to estimate the transformation
window_size = block.InputPort(3).Data;
% length of the data in the previous step
data_length = size(mu_db,1);
if data_length < 5
    % directly use the friction from DB Sensor, if there is insufficient 
    % information to estimate the transformation
    T_FL = [0; 1];
    T_FR = [0; 1];
    T_RL = [0; 1];
    T_RR = [0; 1];
elseif data_length < (window_size + 1)
    % use linear regression to estimate the transformation
    T_FL = [ones( data_length,1 ), mu_db(:,1)]\mu_onboard(:,1);
    T_FR = [ones( data_length,1 ), mu_db(:,2)]\mu_onboard(:,2);
    T_RL = [ones( data_length,1 ), mu_db(:,3)]\mu_onboard(:,3);
    T_RR = [ones( data_length,1 ), mu_db(:,4)]\mu_onboard(:,4);
else
    % use linear regression to estimate the transformation
    T_FL = [ones( window_size,1 ), mu_db( data_length-window_size+1:data_length,1 )]...
        \mu_onboard( data_length-window_size+1:data_length,1 );
    T_FR = [ones( window_size,1 ), mu_db( data_length-window_size+1:data_length,2 )]...
        \mu_onboard( data_length-window_size+1:data_length,2 );
    T_RL = [ones( window_size,1 ), mu_db( data_length-window_size+1:data_length,3 )]...
        \mu_onboard( data_length-window_size+1:data_length,3 );
    T_RR = [ones( window_size,1 ), mu_db( data_length-window_size+1:data_length,4 )]...
        \mu_onboard( data_length-window_size+1:data_length,4 );
end

% maximum data length is limited by the horizon length
if data_length < block.DialogPrm(1).Data
    mu_db = [mu_db; block.InputPort(1).Data'];
    mu_onboard = [mu_onboard; block.InputPort(2).Data'];
else
    mu_db = [mu_db(2:data_length,:); block.InputPort(1).Data'];
    mu_onboard = [mu_onboard(2:data_length,:); block.InputPort(2).Data'];
end

% Transformations to transform DB friction to Onboard friction for all the
% four tyres
transformation = [T_FL'; T_FR'; T_RL'; T_RR'];

block.OutputPort(1).Data = transformation;