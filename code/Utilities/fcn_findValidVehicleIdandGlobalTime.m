%%%%%%%%%%%%%%% Function fcn_findValidVehicleIdandGlobalTime %%%%%%%%%%%%%%
% Purpose:
%   fcn_findValidVehicleIdandGlobalTime queries for unique 
%   (section_id, vehicle_id, global_time) combinations where global_time is 
%   within the bounds defined by timeUnix, timeAimsunLB, and timeAimsunUB. 
%   It also outputs unique (Vehicle ID, Global Time) combinations 
%   representing vehicle trajectories and list of unique section ids.
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   [vehiclePopulation_ID_GlobalTime,SectionId_VehID_GlobalTime,listOfSections] = 
%   fcn_findValidVehicleIdandGlobalTime(timeUnix,timeAimsunLB,timeAimsunUB,names,name_road)
% 
% INPUTS:
%   timeUnix: Unix time
%   timeAimsunLB, timeAimsunUB: Bounds on System Entrance Time
%   names: It's a structure containing name of the database and tables
%   name_road: Name of the road
% 
% OUTPUTS:
%   vehiclePopulation_ID_GlobalTime: It's a Nx2 vector. It consists of
%   unique (vehicle_id, global_time) pairs
%   SectionId_VehID_GlobalTime: It's a Mx3 vector. It consists of unique
%   (section_id, vehicle_id, global_time) pairs
%   listOfSections: It's a Lx1 vector. It consists of unique section_ids
% 
% Examples:
%   timeUnix = 1589072657;
%   timeAimsunLB = 6.61; timeAimsunUB = 1795.61;
%   % name of the road
%   names.road = 'Onekm_straight_highway';
%   % name of the database
%   names.database = 'nsf_roadtraffic_friction_v2';
%   % name of the table that contains the relation between road name and road id
%   names.tableRoad = 'road';
%   % name of the table that contains the relation between road id and section/segment ids
%   names.tableRoadSegment = 'road_road_segment';
%   % name of the table containing traffic data
%   names.tableTraffic = 'road_traffic_procsessed';
%   [vehiclePopulation_ID_GlobalTimeO,SectionId_VehID_GlobalTimeO,listOfSectionsO] = ...
%   fcn_findValidVehicleIdandGlobalTime(timeUnix,timeAimsunLB,timeAimsunUB,names,names.road);
%   [vehiclePopulation_ID_GlobalTimeT,SectionId_VehID_GlobalTimeT,listOfSectionsT] = ...
%   fcn_findValidVehicleIdandGlobalTime(timeUnix,timeAimsunLB,timeAimsunUB,names);
% 
% Author: Satya Prasad
% Created: 2020-05-28
% Updated: 2020-06-23
% 
% Revision history:
% 2020/05/29
%   - Second input is modifed to have all vehicle_id, global_time, and 
%   section_id instead of just section_ids
% 2020/06/02
%   - Modified the inputs from bounds on global time to unix time and 
%   bounds on aimsun time. road_name is replaced with a struct containing
%   road name, database and table names.
% 2020/06/23
%   - Number of inputs are increased to two to have two types of queries. 
%   Added list of unique section ids as one of the output
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vehiclePopulation_ID_GlobalTime,SectionId_VehID_GlobalTime,listOfSections] = ...
    fcn_findValidVehicleIdandGlobalTime(timeUnix,timeAimsunLB,timeAimsunUB,names,name_road)
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
if ~(nargin == 4 || nargin == 5)
    error('fcn_findValidVehicleIdandGlobalTime: Incorrect number of input arguments');
end

%% Query for valid Vehicle ID and Global Time combinations as well as
% Vehicle ID, Global Time, and Section ID combinations
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

% convert table to array
timeUnix     = table2array(timeUnix);
timeAimsunLB = table2array(timeAimsunLB);
timeAimsunUB = table2array(timeAimsunUB);

if nargin == 4
    % SQL query to grab section_id, vehicle_id, and global_time
    sql_query = ['SELECT DISTINCT section_id, vehicle_id, g_time'...
                ' FROM ' names.tableTraffic...
                ' WHERE g_time >= ' num2str(timeUnix+timeAimsunLB)...
                      ' AND g_time <= ' num2str(timeUnix+timeAimsunUB)];
else
    % SQL query to grab section_id, vehicle_id, and global_time
    sql_query = ['SELECT DISTINCT section_id, vehicle_id, g_time'...
                ' FROM ' names.tableTraffic...
                ' WHERE g_time >= ' num2str(timeUnix+timeAimsunLB)...
                      ' AND g_time <= ' num2str(timeUnix+timeAimsunUB)...
                      ' AND section_id IN (SELECT DISTINCT road_segment_id'...
                                         ' FROM ' names.tableRoadSegment...
                                         ' WHERE road_id = (SELECT id'...
                                                          ' FROM ' names.tableRoad...
                                                          ' WHERE name = ''' name_road '''))'];
end
% grab section_id, vehicle_id, and global_time
result = fetch(DB.db_connection, sql_query);
% convert the result from table to array containing
% unique (section_id, vehicle_id, global_time) pairs
SectionId_VehID_GlobalTime = table2array(result);
% unique (vehicle_id, global_time) pairs
vehiclePopulation_ID_GlobalTime = unique(SectionId_VehID_GlobalTime(:,2:3), 'rows');
% unique section_ids
listOfSections = unique(SectionId_VehID_GlobalTime(:,1));

% Disconnect from the database
DB.disconnect();

end