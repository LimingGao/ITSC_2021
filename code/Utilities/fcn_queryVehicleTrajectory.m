%%%%%%%%%%%%%%%%%%%% Function fcn_queryVehicleTrajectory %%%%%%%%%%%%%%%%%%
% Purpose:
%   fcn_queryVehicleTrajectory queries a vehicle's trajectory uniquely 
%   identified by (vehicle_id, global_time) either complete trajectory or
%   just the trajectory on the section 'section_id'
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, ...
%                            global_time, names, section_id)
% 
% INPUTS:
%   vehicle_id: ID of the vehicle
%   global_time: Uniquely identifies a vehicle's trajectory in the database
%   with Vehicle ID
%   names: It's a structure containing name of the database and tables
%   section_id: uniquely identifies a section
% 
% OUTPUTS:
%   queryVehicleTrajectory: It'a a Nx20 matrix. Contains all the 
%   attributes defined by trajectory_attributes sorted in the order of 
%   aimsun_time
% 
% Examples:
%   % Query a valid trajectory
%   vehicle_id = 26; global_time = 1589072708.44;
%   names.database = 'nsf_roadtraffic_friction_v2';
%   names.tableTraffic = 'road_traffic_procsessed';
%   queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, ...
%                            global_time, names);
% 
%   vehicle_id = 26; global_time = 1589072708.44;
%   names.database = 'nsf_roadtraffic_friction_v2';
%   names.tableTraffic = 'road_traffic_procsessed';
%   section_id = 814;
%   queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, ...
%                            global_time, names, section_id);
% 
% Author: Satya Prasad
% Created: 2020-05-27
% Updated: 2020-06-23
% 
% Revision history:
% 2020/05/28
%   - Replaced 'trips_id' with 'g_time'
% 2020/06/02
%   - Added a new input 'names' containing name of the database and
%   name of the table containing traffic data
% 2020/06/18
%   - Refined the input checks. Modified the inputs and query to
%   take section_id into consideration
% ======== to do list ===========
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function queryVehicleTrajectory = fcn_queryVehicleTrajectory(vehicle_id, ...
                                  global_time, names, section_id)
trajectory_attributes = ['vehicle_id, vehicle_type, vehicle_length, vehicle_width, '...
                        'section_id, junction_id, lane_number, direction, '...
                        'latitude, longitude, altitude, '...
                        'positioncg_e, positioncg_n, positioncg_u, yaw, current_speed, '...
                        'station_section, station_total, aimsun_time, g_time']; % attributes in the vehicle trajectory

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
if ~(nargin == 3 || nargin == 4)
    error('fcn_queryVehicleTrajectory: Incorrect number of input arguments.')
end

% Check the size and validity of inputs
if ~isnumeric(vehicle_id) || (length(vehicle_id) ~= 1) || ...
        any(vehicle_id <= 0) || any(vehicle_id ~= round(vehicle_id))
    % display an error message if 'vehicle_id' is not a positive integer
    error('vehicle_id must be a POSITIVE INTEGER')
    
elseif ~isnumeric(global_time) || (length(global_time) ~= 1) || ...
        any(global_time <= 0)
    % display an error message if 'global_time' is not positive
    error('global_time must be POSITIVE')
    
end

%% Query vehicle trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create an instance
% connect to the database
DB = Database(names.database);

if nargin == 4
    % Check the size and validity of inputs
    if ~isnumeric(section_id) || (length(section_id) ~= 1) || ...
            any(section_id <= 0) || any(section_id ~= round(section_id))
        % display an error message if 'section_id' is not a positive integer
        error('section_id must be a POSITIVE INTEGER')
    end
    % SQL query to grab vehicle trajectory in a section
    traj_query = ['SELECT ' trajectory_attributes...
                 ' FROM ' names.tableTraffic...
                 ' WHERE vehicle_id = ' num2str(vehicle_id) ...
                       ' AND g_time = ' num2str(global_time) ...
                       ' AND section_id = ' num2str(section_id) ...
                 ' ORDER BY aimsun_time'];
    
else
    % SQL query to grab vehicle trajectory
    traj_query = ['SELECT ' trajectory_attributes...
                 ' FROM ' names.tableTraffic...
                 ' WHERE vehicle_id = ' num2str(vehicle_id) ...
                       ' AND g_time = ' num2str(global_time) ...
                 ' ORDER BY aimsun_time'];
end

% grab trajectory data from the DB
result = fetch(DB.db_connection, traj_query);
% convert the result from table to array
queryVehicleTrajectory = table2array(result);

% Disconnect from the database
DB.disconnect();

end