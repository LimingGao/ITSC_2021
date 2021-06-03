%%%%%%%%%%%%%%% S-Function sfcn_dbToOnboardTransformation %%%%%%%%%%%%%%%%
% Purpose:
%   Estimates the transformation to transform friction value from DB sensor
%   to a value from Onboard sensor
% 
% INPUTS:
%   INPUT 1: friction value from the DB sensor
%   Input 2: friction value from the Onboard sensor
% 
% OUTPUTS:
%   OUTPUT 1: intercept and slope of the transformation estimate
% 
% PARAMETERS:
%   Parameter 1: window size to estimate transformation
% 
% Author: Satya Prasad
% Create Date: 2020-10-01
% Revision history:
% 2020\10\02: The range used to estimate the transformation is limited to
% the last step.
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_dbToOnboardTransformation(block)

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

% Register the parameters
block.NumDialogPrms     = 1;
block.DialogPrmsTunable = {'Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = 4;
block.InputPort(2).Dimensions = 4;

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
persistent onboard_mu_FL onboard_mu_FR onboard_mu_RL onboard_mu_RR;
% variables to store friction values from the DB sensor
persistent db_mu_FL db_mu_FR db_mu_RL db_mu_RR;
% variables to store intercept and slope of transformation
persistent transformation

if block.Dwork(1).Data
    % Intialize the variables
    db_mu_FL = [];
    db_mu_FR = [];
    db_mu_RL = [];
    db_mu_RR = [];
    onboard_mu_FL = [];
    onboard_mu_FR = [];
    onboard_mu_RL = [];
    onboard_mu_RR = [];
    
    % reset the variable
    block.Dwork(1).Data = 0;
    
end

if length(db_mu_FL) < 5
    % directly use the friction from DB Sensor, if there is insufficient 
    % information to estimate the transformation
    db_mu_FL = [db_mu_FL; block.InputPort(1).Data(1)];
    db_mu_FR = [db_mu_FR; block.InputPort(1).Data(2)];
    db_mu_RL = [db_mu_RL; block.InputPort(1).Data(3)];
    db_mu_RR = [db_mu_RR; block.InputPort(1).Data(4)];
    onboard_mu_FL = [onboard_mu_FL; block.InputPort(2).Data(1)];
    onboard_mu_FR = [onboard_mu_FR; block.InputPort(2).Data(2)];
    onboard_mu_RL = [onboard_mu_RL; block.InputPort(2).Data(3)];
    onboard_mu_RR = [onboard_mu_RR; block.InputPort(2).Data(4)];
    
    T_FL = [0; 1];
    T_FR = [0; 1];
    T_RL = [0; 1];
    T_RR = [0; 1];
    
else
    % use linear regression to estimate the transformation
    T_FL = [ones(size(db_mu_FL)), db_mu_FL]\onboard_mu_FL;
    T_FR = [ones(size(db_mu_FR)), db_mu_FR]\onboard_mu_FR;
    T_RL = [ones(size(db_mu_RL)), db_mu_RL]\onboard_mu_RL;
    T_RR = [ones(size(db_mu_RR)), db_mu_RR]\onboard_mu_RR;
    
    % maximum size of the data is limited by the first parameter
    if length(db_mu_FL) < block.DialogPrm(1).Data
        db_mu_FL = [db_mu_FL; block.InputPort(1).Data(1)];
        db_mu_FR = [db_mu_FR; block.InputPort(1).Data(2)];
        db_mu_RL = [db_mu_RL; block.InputPort(1).Data(3)];
        db_mu_RR = [db_mu_RR; block.InputPort(1).Data(4)];
        onboard_mu_FL = [onboard_mu_FL; block.InputPort(2).Data(1)];
        onboard_mu_FR = [onboard_mu_FR; block.InputPort(2).Data(2)];
        onboard_mu_RL = [onboard_mu_RL; block.InputPort(2).Data(3)];
        onboard_mu_RR = [onboard_mu_RR; block.InputPort(2).Data(4)];
        
    else
        db_mu_FL = [db_mu_FL(2:end); block.InputPort(1).Data(1)];
        db_mu_FR = [db_mu_FR(2:end); block.InputPort(1).Data(2)];
        db_mu_RL = [db_mu_RL(2:end); block.InputPort(1).Data(3)];
        db_mu_RR = [db_mu_RR(2:end); block.InputPort(1).Data(4)];
        onboard_mu_FL = [onboard_mu_FL(2:end); block.InputPort(2).Data(1)];
        onboard_mu_FR = [onboard_mu_FR(2:end); block.InputPort(2).Data(2)];
        onboard_mu_RL = [onboard_mu_RL(2:end); block.InputPort(2).Data(3)];
        onboard_mu_RR = [onboard_mu_RR(2:end); block.InputPort(2).Data(4)];
        
    end
    
end

% Transformations to transform DB friction to Onboard friction for all the
% four tyres
transformation(1,:) = T_FL;
transformation(2,:) = T_FR;
transformation(3,:) = T_RL;
transformation(4,:) = T_RR;

block.OutputPort(1).Data = transformation;