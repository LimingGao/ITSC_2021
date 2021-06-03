%%%%%%%%%%%%%% Function fcn_swervedVehicleTrajectory %%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   fcn_swervedVehicleTrajectory adds sine wave distortion to 
%   queryVehicleTrajectory iff flag_swerve is true
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   swervedVehicleTrajectory = fcn_swervedVehicleTrajectory(queryVehicleTrajectory,...
%                      amplitude, time_period, refLLA, fields, flag_swerve)
% 
% INPUTS:
%   queryVehicleTrajectory: Nx20 matrix where N is atleast 1.
% 
%   amplitude: It is either a row or column vector.
% 
%   time_period: It is either a row or column vector. amplitude and
%   time_period must have same length.
% 
%   lateral_offset: It is a scalar to add a fixed lateral offset to the
%   vehicle trajectory.
% 
%   refLLA: reference LLA for ENU to LLA transformation
% 
%   fields: structure mapping column numbers of queryVehicleTrajectory to column names
% 
%   flag_swerve: 'true' indicates add distortion. 'false' indicates no
%   distortion
% 
% OUTPUTS:
%   swervedVehicleTrajectory: Contains swerved trajectory
% 
% Author: Satya Prasad
% Create Date: 2020-05-28
% Revision history:
% 2020/06/02: Added reference LLA as an input and also structure containing
% column numbers as an input
% 2020/08/03: Added fixed 'lateral_offset' as one of the input. Modified
% the if statements to make the code compact.
% 
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function swervedVehicleTrajectory = fcn_swervedVehicleTrajectory(vehicleTrajectory,...
                       amplitude, time_period, lateral_offset, refLLA, fields, flag_swerve)
do_debug = 0; % DEBUG flag. set it to one to compare input and output trajectories

%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are there right number of inputs?
switch nargin
    case 7
    
    case 6
        flag_swerve = true;
    otherwise
        error('Incorrect number of input arguments');
end

% Check the size and validity of inputs
if (length(amplitude) ~= length(time_period))
    % display an error message if 'amplitude' and 'time_period' are of 
    % different length
    error('ERROR: amplitude and time_period must be of same size');
    
elseif size(vehicleTrajectory,2) ~= 20
    % display an error message if 'queryVehicleTrajectory' has different
    % number of columns than expected
    error('ERROR: queryVehicleTrajectory must be a Nx20 matrix');
    
elseif ~isa(flag_swerve, 'logical')
    error('ERROR: flag_swerve must be a bool');
end

%% Add swerving to vehicleTrajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, unique_indices] = unique(vehicleTrajectory(:, fields.station)); % Indices of unique station coordinates
vehicleTrajectory   = vehicleTrajectory(unique_indices, :); % Trajectory with unique station coordinates

path = vehicleTrajectory(:, [fields.east, fields.north]); % Path (East and North coordinates) of the trajectory
vehicleTrajectory(:, fields.station) = [0; cumsum(sqrt(sum(diff(path).^2, 2)))]; % Update the station coordinates based on location of the vehicle

num_points = size(vehicleTrajectory, 1); % the default number of points to use

%% Fill in the array of reference station points
reference_traversal = fcn_Path_convertPathToTraversalStructure(path);

%% Fill in the array of offset distances.
offset_from_reference = zeros(num_points, 1); % Initialization
if flag_swerve
    swerving_time = 0.1*(1:num_points)'; % time scale used to add sine distortion
    for i = 1:length(amplitude)
        offset_from_reference = offset_from_reference+...
            amplitude(i)*sin(2*pi*swerving_time/time_period(i));
    end
    offset_from_reference = offset_from_reference+lateral_offset;
else
    offset_from_reference = offset_from_reference+lateral_offset; % Only lateral-offset to the trajectory
end % NOTE: Ends flag_swerve

%% Find projection from reference orthogonally
% Set the projection type (see help in function below for details)
flag_rounding_type = 4; % This averages the projection vectors along segments

% Find the unit normal vectors at each of the station points
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    reference_traversal.Station, reference_traversal, flag_rounding_type);
unit_vectors = unit_normal_vector_end-unit_normal_vector_start;

%% Perform iterations over trajectories
offset_path = unit_normal_vector_start + unit_vectors.*offset_from_reference;

%% Update the output trajectory
% Position in ENU and LLA, Station and Yaw
offset_east    = offset_path(1);
offset_north   = offset_path(2);
offset_up      = vehicleTrajectory(:, fields.up);
offset_station = [0; cumsum(sqrt(sum(diff(offset_path).^2, 2)))];
offset_yaw     = fcn_Path_calcYawFromPathSegments(offset_path);
[offset_lat, offset_lon, offset_alt] = ...
    enu2geodetic(offset_east, offset_north, offset_up, ...
    refLLA(1), refLLA(2), refLLA(3), referenceEllipsoid('wgs84')); % convert from ENU to LLA

swervedVehicleTrajectory = vehicleTrajectory; % Intialize the ouput trajectory with input trajectory of unique station coordinates
swervedVehicleTrajectory(:, fields.east)    = offset_east;
swervedVehicleTrajectory(:, fields.north)   = offset_north;
swervedVehicleTrajectory(:, fields.station) = offset_station;
swervedVehicleTrajectory(:, fields.yaw)     = [offset_yaw; offset_yaw(end)];
swervedVehicleTrajectory(:, fields.lat)     = offset_lat;
swervedVehicleTrajectory(:, fields.lon)     = offset_lon;
swervedVehicleTrajectory(:, fields.alt)     = offset_alt;

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_debug
    figure(92345)
    clf
    plot(vehicleTrajectory(:, fields.east), vehicleTrajectory(:, fields.north), 'b*');
    hold on;
    plot(swervedVehicleTrajectory(:, fields.east), ...
        swervedVehicleTrajectory(:, fields.north), 'gd');
    title('Vehicle Trajectory');
    xlabel('east(m)'); ylabel('north(m)');
    legend('Reference Trajectory', 'Swerved Trajectory');
end

end