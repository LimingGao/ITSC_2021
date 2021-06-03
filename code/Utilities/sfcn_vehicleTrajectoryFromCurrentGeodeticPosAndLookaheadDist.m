%%%%%%%%% S-Function sfcn_vehicleTrajectoryFromCurrentGeodeticPosAndLookaheadDist %%%%%%%%%
% Purpose:
%   Finds the point on the interpolated trajectory that is at a
%   look_ahead_distance from the current vehicle position
% 
% Distance metric is 'Spherical Law of Cosines'
% 
% INPUTS:
%   INPUT 1: 'input_trajectory' is the complete vehicle_trajectory
%   available in the database
%   INPUT 2: 'vehicle_cg' is the cg of the vehicle in LLA
%   INPUT 3: look ahead distance
% 
% OUTPUTS:
%   OUTPUT 1: position of the vehicle ahead of Input 1 at a
%   look_head_distance
% 
% PARAMETERS:
%   Parameter 1: Size of INPUT 1
%   Parameter 2: Parameter limiting the search range for nearest neighbour
%   Parameter 3: Structure mapping column names of INPUT 1 with column
%   numbers
% 
% Author: Satya Prasad
% Created: 2020-06-25
% 
% Revision History:
% 2020/06/27
%   - Commented out the distance computation using spherical law of
%   cosines. Added the haversine formula to compute distance.
%   'flag_full_window' is made as a parameter.
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sfcn_vehicleTrajectoryFromCurrentGeodeticPosAndLookaheadDist(block)
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 3;
block.NumOutputPorts = 1;
% Register the parameters
block.NumDialogPrms     = 4;
block.DialogPrmsTunable = {'Nontunable', 'Nontunable', 'Nontunable', 'Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = block.DialogPrm(1).Data;
block.InputPort(2).Dimensions = 3;
block.InputPort(3).Dimensions = 1;

% Override output port properties
block.OutputPort(1).Dimensions = block.DialogPrm(1).Data(2);

% Set block sample time to be inherited
block.SampleTimes = [-1, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('Outputs', @Output);   % Required

function Output(block)
%% Output function with search window depending on vehicle trajectory
% 
% 'lad' - look_ahead_distance
% 
% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
% variable to store trajectory with unique station coordinates
persistent input_trajectory;
% size of input trajectory
persistent size_input_trajectory;
% variable to do the initialization
persistent flag;
% variable to store sine and cosine of latitude of points on the trajectory
persistent sineAndCosineOfLat
% variable to store lower bounds and window size to limit the search range 
% for nearest neighbour
persistent start_index start_index_lad half_window_size;

% current cg of the vehicle in LLA
vehicle_cg = block.InputPort(2).Data;
cg_lat = vehicle_cg(1);
cg_lon = vehicle_cg(2);
cg_alt = vehicle_cg(3);

refLat = -48.876667; refLon = -123.393333; refAlt = 0;  % reference LLA
R = 6371e3; % earth's radius in meters
% look ahead distance
look_ahead_distance = block.InputPort(3).Data;
% structure mapping column names to column numbers
fieldsPreSimulationVehicleTrajectory = block.DialogPrm(3).Data;
% flag to switch between only looking forward and looking in both the
% directions; 'true' -> both direction; 'false' -> only forward
flag_full_window = block.DialogPrm(4).Data;

% set the flag to empty if there is a new trajectory
if ~isempty(input_trajectory)
    if ~any(input_trajectory(1,[fieldsPreSimulationVehicleTrajectory.vehId,...
            fieldsPreSimulationVehicleTrajectory.globalTime]) == ...
            block.InputPort(1).Data(1,[fieldsPreSimulationVehicleTrajectory.vehId,...
            fieldsPreSimulationVehicleTrajectory.globalTime]))
        flag = [];
    end
end

if isempty(flag)
    %% Initialization
    % find all the indices, on the trajectory, that have unique station coordinates
    [~,unique_indices] = unique(block.InputPort(1).Data(:,...
                            fieldsPreSimulationVehicleTrajectory.station));
    % sort the unique_indices in the increasing order
    sorted_indices     = sortrows(unique_indices);
    % grab the trajectory that have unique station coordinates
    input_trajectory   = block.InputPort(1).Data(sorted_indices,:);
    size_input_trajectory = size(input_trajectory);
    % sine and cosine of latitude of points on the trajectory
    sineAndCosineOfLat = [sin(input_trajectory(:,fieldsPreSimulationVehicleTrajectory.lat)),...
                          cos(input_trajectory(:,fieldsPreSimulationVehicleTrajectory.lat))];
    % variable to limit the search range
    half_window_size = max(2,ceil(block.DialogPrm(2).Data/min(diff(...
      input_trajectory(:,fieldsPreSimulationVehicleTrajectory.station)))));
    
    
    %% Nearest neighbour search based on current vehicle position
    % current vehicle position is an input to the function
    
%     % compute the distance using spherical law of cosines
%     cosinelaw_distance_between_trajectory_and_cg = R*acos(sin(cg_lat)*sineAndCosineOfLat(:,1) + ...
%                                                           cos(cg_lat)*sineAndCosineOfLat(:,2).*...
%                     cos(input_trajectory(:,fieldsPreSimulationVehicleTrajectory.lon)-cg_lon));
%     % pick the indices of smallest two elements
%     [~,nearest_index] = mink(cosinelaw_distance_between_trajectory_and_cg,2);
    % compute the distance using haversine law of cosines
    haversine_a = sin(0.5*(input_trajectory(:,fieldsPreSimulationVehicleTrajectory.lat)-cg_lat)).^2 + ...
                  cos(cg_lat)*sineAndCosineOfLat(:,2).*...
                  sin(0.5*(input_trajectory(:,fieldsPreSimulationVehicleTrajectory.lon)-cg_lon)).^2;
    haversine_c = 2*atan2(sqrt(haversine_a), sqrt(1-haversine_a));
    haversine_distance_between_trajectory_and_cg = R*haversine_c;
    % pick the indices of smallest two elements
    [~, nearest_index] = mink(haversine_distance_between_trajectory_and_cg,2);
    % sort the indices in the increasing order
    nearest_index     = sortrows(nearest_index);
    
    % intialize the lower bound of the search window
    if flag_full_window
        start_index = max(1,nearest_index(1)-half_window_size);
    else
        start_index = nearest_index(1);
    end
    
    if size(nearest_index,1) >= 2
        % index of nearest point, on trajectory, behind the input position
        rear_index  = nearest_index(1);
        % index of nearest point, on trajectory, infront of the input position
        front_index = nearest_index(2);
        
        % nearest point, on trajectory, behind the vehicle_cg
        lat_rear     = input_trajectory(rear_index, fieldsPreSimulationVehicleTrajectory.lat);
        lon_rear     = input_trajectory(rear_index, fieldsPreSimulationVehicleTrajectory.lon);
        alt_rear     = input_trajectory(rear_index, fieldsPreSimulationVehicleTrajectory.alt);
        station_rear = input_trajectory(rear_index, fieldsPreSimulationVehicleTrajectory.station);
        
        % nearest point, on trajectory, infront of the vehicle_cg
        lat_front     = input_trajectory(front_index, fieldsPreSimulationVehicleTrajectory.lat);
        lon_front     = input_trajectory(front_index, fieldsPreSimulationVehicleTrajectory.lon);
        alt_front     = input_trajectory(front_index, fieldsPreSimulationVehicleTrajectory.alt);
        station_front = input_trajectory(front_index, fieldsPreSimulationVehicleTrajectory.station);
        
    elseif size(nearest_index,1) == 1
        error('ERROR: Only one nearest neighbour is available \n');
    else
        error('ERROR: There is no nearest neighbour \n');
    end
    
    arrayLat = [cg_lat; lat_front];
    arrayLon = [cg_lon; lon_front];
    arrayAlt = [cg_alt; alt_front];
    enu = fcn_LLA2ENU(arrayLat, arrayLon, arrayAlt, lat_rear, lon_rear, alt_rear);
    % vector pointing from trajectory towards the vehicle
    point_vector = enu(1,1:2);
    % vector pointing, along the trajectory, in the direction of travel
    path_vector  = enu(2,1:2);
    
    % station of the nearest neigbour at a look_ahead_distance, on the
    % interpolated trajectory, to the vehicle_cg
    multiplying_factor = dot(path_vector, point_vector)/(norm(path_vector)^2);
    
    cg_station_lad = station_rear ...
                     + (station_front-station_rear)*multiplying_factor ...
                     + look_ahead_distance;
    
    %% Nearest neighbour search based on cg_station_lad
    
    % compute absolute difference between station coordinates of every
    % point in input_trajectory and cg_station_lad
    diff_between_trajectory_and_cg_lad_stations = abs(input_trajectory(:,fieldsPreSimulationVehicleTrajectory.station) ...
                                                  - cg_station_lad);
    % pick the indices of smallest two elements
    [~,nearest_index_lad] = mink(diff_between_trajectory_and_cg_lad_stations,2);
    % sort the indices in the increasing order
    nearest_index_lad     = sortrows(nearest_index_lad);
    
    % index that decides ouput properties other than position and station
    output_index = nearest_index_lad(1);
    % initialize the lower bound for search window
    if flag_full_window
        start_index_lad = max(1,nearest_index_lad(1)-half_window_size);
    else
        start_index_lad = nearest_index_lad(1);
    end
    
    if size(nearest_index_lad,1) >= 2
        % index of nearest point, on trajectory, behind the desired station
        rear_index_lad  = nearest_index_lad(1);
        % index of nearest point, on trajectory, infront of the desired station
        front_index_lad = nearest_index_lad(2);
        
        % nearest point, on trajectory, behind the desired station
        lat_rear_lad     = input_trajectory(rear_index_lad, fieldsPreSimulationVehicleTrajectory.lat);
        lon_rear_lad     = input_trajectory(rear_index_lad, fieldsPreSimulationVehicleTrajectory.lon);
        alt_rear_lad     = input_trajectory(rear_index_lad, fieldsPreSimulationVehicleTrajectory.alt);
        station_rear_lad = input_trajectory(rear_index_lad, fieldsPreSimulationVehicleTrajectory.station);
        
        % nearest point, on trajectory, infront of the desired station
        lat_front_lad     = input_trajectory(front_index_lad, fieldsPreSimulationVehicleTrajectory.lat);
        lon_front_lad     = input_trajectory(front_index_lad, fieldsPreSimulationVehicleTrajectory.lon);
        alt_front_lad     = input_trajectory(front_index_lad, fieldsPreSimulationVehicleTrajectory.alt);
        station_front_lad = input_trajectory(front_index_lad, fieldsPreSimulationVehicleTrajectory.station);
        
    elseif size(nearest_index_lad,1) == 1
        error('ERROR: Only one nearest neighbour is available \n');
    else
        error('ERROR: There is no nearest neighbour \n');
    end
    
    % vector pointing, along the trajectory, in the direction of travel
    path_vector_lad = [lat_front_lad-lat_rear_lad, lon_front_lad-lon_rear_lad, alt_front_lad-alt_rear_lad];
    
    % position of the nearest neigbour, on the interpolated trajectory, to
    % the cg_station_lad
    multiplying_factor_lad = (cg_station_lad-station_rear_lad)/(station_front_lad-station_rear_lad);
    
    cg_lat_lad = lat_rear_lad + path_vector_lad(1)*multiplying_factor_lad;
    cg_lon_lad = lon_rear_lad + path_vector_lad(2)*multiplying_factor_lad;
    cg_alt_lad = alt_rear_lad + path_vector_lad(3)*multiplying_factor_lad;
    
    % convert cg of the vehicle to ENU
    enu_lad = fcn_LLA2ENU(cg_lat_lad,cg_lon_lad,cg_alt_lad,refLat,refLon,refAlt);
    cg_east_lad  = enu_lad(1);
    cg_north_lad = enu_lad(2);
    cg_up_lad    = enu_lad(3);
    
    flag = 1;
    
else
    %% Nearest neighbour search based on current vehicle position
    % current vehicle position is an input to the function
    
    % if start_index reaches maximum, then set it to one less than the maximum
    if start_index == size_input_trajectory(1)
        start_index = size_input_trajectory(1) - 1;
    end
    
    % set the upper bound to limit the search range for nearest neighbour
    % depending on size of INPUT 1/Parameter 1 and half_window_size
    if flag_full_window
        end_index = min(start_index+2*half_window_size+1,size_input_trajectory(1));
    else
        end_index = min(start_index+half_window_size,size_input_trajectory(1));
    end
    
    % select the trajectory that lies within bounds
    selected_trajectory = input_trajectory(start_index:end_index, ...
                        [fieldsPreSimulationVehicleTrajectory.lat, ...
                         fieldsPreSimulationVehicleTrajectory.lon, ...
                         fieldsPreSimulationVehicleTrajectory.alt, ...
                         fieldsPreSimulationVehicleTrajectory.station]);
    % select the sine and cosine of latitude of points on the trajectory that lies within bounds
    selected_sineAndCosineOfLat = sineAndCosineOfLat(start_index:end_index,:);
    
%     % compute the distance using spherical law of cosines
%     cosinelaw_distance_between_trajectory_and_cg = R*acos(sin(cg_lat)*selected_sineAndCosineOfLat(:,1) + ...
%                                                           cos(cg_lat)*selected_sineAndCosineOfLat(:,2).*...
%                                                           cos(selected_trajectory(:,2)-cg_lon));
%     % pick the indices of smallest two elements
%     [~, nearest_index] = mink(cosinelaw_distance_between_trajectory_and_cg,2);
    % compute the distance using haversine law of cosines
    haversine_a = sin(0.5*(selected_trajectory(:,1)-cg_lat)).^2 + ...
                  cos(cg_lat)*selected_sineAndCosineOfLat(:,2).*...
                  sin(0.5*(selected_trajectory(:,2)-cg_lon)).^2;
    haversine_c = 2*atan2(sqrt(haversine_a), sqrt(1-haversine_a));
    haversine_distance_between_trajectory_and_cg = R*haversine_c;
    % pick the indices of smallest two elements
    [~, nearest_index] = mink(haversine_distance_between_trajectory_and_cg,2);
    % sort the indices in the increasing order
    nearest_index      = sortrows(nearest_index);
    
    % shift start_index depending on the nearest neighbour
    if flag_full_window
        start_index = max(1,start_index+nearest_index(1)-1-half_window_size);
    else
        start_index = start_index+nearest_index(1)-1;
    end
    
    if size(nearest_index,1) >= 2
        % index of nearest point, on trajectory, behind the input position
        rear_index  = nearest_index(1);
        % index of nearest point, on trajectory, infront of the input position
        front_index = nearest_index(2);
        
        % nearest point, on trajectory, behind the vehicle_cg
        lat_rear     = selected_trajectory(rear_index, 1);
        lon_rear     = selected_trajectory(rear_index, 2);
        alt_rear     = selected_trajectory(rear_index, 3);
        station_rear = selected_trajectory(rear_index, 4);
        
        % nearest point, on trajectory, infront of the vehicle_cg
        lat_front     = selected_trajectory(front_index, 1);
        lon_front     = selected_trajectory(front_index, 2);
        alt_front     = selected_trajectory(front_index, 3);
        station_front = selected_trajectory(front_index, 4);
        
    elseif size(nearest_index,1) == 1
        error('ERROR: Only one nearest neighbour is available \n');
    else
        error('ERROR: There is no nearest neighbour \n');
    end
    
    arrayLat = [cg_lat; lat_front];
    arrayLon = [cg_lon; lon_front];
    arrayAlt = [cg_alt; alt_front];
    enu = fcn_LLA2ENU(arrayLat, arrayLon, arrayAlt, lat_rear, lon_rear, alt_rear);
    % vector pointing from trajectory towards the vehicle
    point_vector = enu(1,1:2);
    % vector pointing, along the trajectory, in the direction of travel
    path_vector  = enu(2,1:2);
    
    % station of the nearest neigbour at a look_ahead_distance, on the
    % interpolated trajectory, to the vehicle_cg
    multiplying_factor = dot(path_vector, point_vector)/(norm(path_vector)^2);
    
    cg_station_lad = station_rear ...
                     + (station_front-station_rear)*multiplying_factor ...
                     + look_ahead_distance;
    
    %% Nearest neighbour search based on cg_station_lad
    
    % if start_index_lad reaches maximum, then set it to one less than the maximum
    if start_index_lad == size_input_trajectory(1)
        start_index_lad = size_input_trajectory(1) - 1;
    end
    
    % set the upper bound to limit the search range for nearest neighbour
    % depending on size of INPUT 1/Parameter 1 and half_window_size
    if flag_full_window
        end_index_lad = min(start_index_lad+2*half_window_size+1,size_input_trajectory(1));
    else
        end_index_lad = min(start_index_lad+half_window_size,size_input_trajectory(1));
    end
    
    % select the trajectory that lies within bounds
    selected_trajectory_lad = input_trajectory(start_index_lad:end_index_lad, ...
                            [fieldsPreSimulationVehicleTrajectory.lat, ...
                             fieldsPreSimulationVehicleTrajectory.lon, ...
                             fieldsPreSimulationVehicleTrajectory.alt, ...
                             fieldsPreSimulationVehicleTrajectory.station]);
    
    % compute absolute difference between station coordinates of every point
    % in selected_trajectory_lad and cg_station_lad
    diff_between_trajectory_and_cg_lad_stations = abs(selected_trajectory_lad(:,4) ...
                                                  - cg_station_lad);
    % pick the indices of smallest two elements
    [~, nearest_index_lad] = mink(diff_between_trajectory_and_cg_lad_stations, 2);
    % sort the indices in the increasing order
    nearest_index_lad      = sortrows(nearest_index_lad);
    
    % index that decides ouput properties other than position and station
    output_index = start_index_lad + nearest_index_lad(1) - 1;
    % shift start_index_lad depending on the nearest neighbour
    if flag_full_window
        start_index_lad = max(1,start_index_lad+nearest_index_lad(1)-1-half_window_size);
    else
        start_index_lad = start_index_lad+nearest_index_lad(1)-1;
    end
    
    if size(nearest_index_lad,1) >= 2
        % index of nearest point, on trajectory, behind the desired station
        rear_index_lad  = nearest_index_lad(1);
        % index of nearest point, on trajectory, infront of the desired station
        front_index_lad = nearest_index_lad(2);
        
        % nearest point, on trajectory, behind the desired station
        lat_rear_lad     = selected_trajectory_lad(rear_index_lad, 1);
        lon_rear_lad     = selected_trajectory_lad(rear_index_lad, 2);
        alt_rear_lad     = selected_trajectory_lad(rear_index_lad, 3);
        station_rear_lad = selected_trajectory_lad(rear_index_lad, 4);
        
        % nearest point, on trajectory, infront of the desired station
        lat_front_lad     = selected_trajectory_lad(front_index_lad, 1);
        lon_front_lad     = selected_trajectory_lad(front_index_lad, 2);
        alt_front_lad     = selected_trajectory_lad(front_index_lad, 3);
        station_front_lad = selected_trajectory_lad(front_index_lad, 4);
        
    elseif size(nearest_index_lad,1) == 1
        error('ERROR: Only one nearest neighbour is available \n');
    else
        error('ERROR: There is no nearest neighbour \n');
    end
    
    % vector pointing, along the trajectory, in the direction of travel
    path_vector_lad  = [lat_front_lad-lat_rear_lad, lon_front_lad-lon_rear_lad, alt_front_lad-alt_rear_lad];
    
    % position of the nearest neigbour, on the interpolated trajectory, to
    % the cg_station_lad
    multiplying_factor_lad = (cg_station_lad-station_rear_lad)/(station_front_lad-station_rear_lad);
    
    cg_lat_lad = lat_rear_lad + path_vector_lad(1)*multiplying_factor_lad;
    cg_lon_lad = lon_rear_lad + path_vector_lad(2)*multiplying_factor_lad;
    cg_alt_lad = alt_rear_lad + path_vector_lad(3)*multiplying_factor_lad;
    % convert the cg to ENU
    enu_lad = fcn_LLA2ENU(cg_lat_lad,cg_lon_lad,cg_alt_lad,refLat,refLon,refAlt);
    cg_east_lad  = enu_lad(1);
    cg_north_lad = enu_lad(2);
    cg_up_lad    = enu_lad(3);
    
end

% initialize the output value
vehicle_trajectory = input_trajectory(output_index,:);
% update the output with estimated values
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.east)    = cg_east_lad;
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.north)   = cg_north_lad;
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.up)      = cg_up_lad;
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.station) = cg_station_lad;
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.lat)     = cg_lat_lad;
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.lon)     = cg_lon_lad;
vehicle_trajectory(fieldsPreSimulationVehicleTrajectory.up)      = cg_alt_lad;

block.OutputPort(1).Data = vehicle_trajectory;