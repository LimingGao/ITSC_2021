%%%%%%%%%%%%%%%%%%% S-Function sfcn_windowSizeUsingAVAR %%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Estimates the window size for computing Moving Average using AVAR.
% 
% INPUTS:
%   Input 1: friction value from the Onboard sensor.
% 
% OUTPUTS:
%   OUTPUT 1: window size for computing Moving Average.
% 
% PARAMETERS:
%   PARAMETER 1: horizon length - decides the maximum length of the data to
%   be stored.
% 
% Author: Satya Prasad
% Create Date: 2020-10-06
% Revision history:
% 2020\10\06: 
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_windowSizeUsingAVAR(block)

% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

function setup(block)
%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
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
block.InputPort(1).Dimensions = 4;

% Override output port properties
block.OutputPort(1).Dimensions = 1;

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
%% Output function to estimate window size using AVAR

% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
persistent mu

if block.Dwork(1).Data
    % Intialize the variables for every new simulation run
    mu = [];
    % reset the variable
    block.Dwork(1).Data = 0;
end

% total number of sensor measurements
M = numel( mu );
if M == 0
    p = -1;
else
    % find the power of 2 that is closest as well as less than the total 
    % number of measurements
    p = floor( log2( M ) ); 
end

if p >= 4
    % number of elements after truncation
    n = 2^p;
    % truncate measurements to the nearest power of 2
    mu_truncated = mu( (M-n+1):M );
    % variable to store allan variance
    avar = zeros( min( (p-2), 9 ), 1 );
    % variable to store data block size
    tau  = zeros( min( (p-2), 9 ), 1 );
    for i = 2:min( (p-1), 10 )
        m = 2^i; % data block size
        avar_sum = 0;
        for k = (2*m):n
            avar_sum = avar_sum + (mean( mu_truncated( (k-m+1):k ) )...
                                 - mean( mu_truncated( (k-2*m+1):(k-m) ) ))^2;
        end
        avar(i-1) = 0.5*avar_sum/(n - (2*m)); % write avar to the output
        tau(i-1)  = m; % write data block size to the output
    end
    [~,ind] = min( avar ); % find the index corresponding to minimum avar
    window_size  = tau(ind); % block size corresponding to minimum avar
else
    window_size = 4; % default block size based on avar is set to '4'
end

% Add new data after estimating the window size
if numel( mu ) < block.DialogPrm(1).Data
    mu = [mu; block.InputPort(1).Data(1)];
else
    mu = [mu(2:end); block.InputPort(1).Data(1)];
end

block.OutputPort(1).Data = window_size;