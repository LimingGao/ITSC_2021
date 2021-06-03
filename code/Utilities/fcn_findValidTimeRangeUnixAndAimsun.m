%%%%%%%%%%%%%%% Function fcn_findValidTimeRangeUnixAndAimsun %%%%%%%%%%%%%%
% Purpose:
%   fcn_findValidTimeRangeUnixAndAimsun queries valid trip id, trip description, 
%   unix time, and bounds on system entrance time related to either all 
%   trips or trips on 'name_road'
% 
% Matlab work Path: ~\GitHub\forgetfulDBs\Utilities
% 
% Format:
%   unixTimeAndAimsunTimeRange = fcn_findValidTimeRangeUnixAndAimsun(names,name_road)
% 
% INPUTS:
%   names: It is a structure containing names of database and relations.
%   name_road: Name of the road.
% 
% OUTPUTS:
%   unixTimeAndAimsunTimeRange: It's a Nx5 table. First and second columns
%   contian trip id and trip description. Third column contains Unix Time. 
%   Fourth and fifth columns contains lower and upper bounds on System 
%   Entrance Time respectively. Output contains information of all trips if 
%   number of inputs is one. Output contains information of all trips on 
%   'name_road' if number of inputs is two.
% 
% Examples:
%   % name of the road
%   names.road = 'Onekm_straight_highway';
%   % name of the database
%   names.database = 'nsf_roadtraffic_friction_v2';
%   % name of the table containing trip information
%   names.tableTrips = 'trips';
%   % name of the table that contains the relation between road name and road id
%   names.tableRoad = 'road';
%   % name of the table that contains the relation between road id and section/segment ids
%   names.tableRoadSegment = 'road_road_segment';
%   % name of the table containing traffic data
%   names.tableTraffic = 'road_traffic_procsessed';
%   unixTimeAndAimsunTimeRangeO = fcn_findValidTimeRangeUnixAndAimsun(names,names.road);
%   unixTimeAndAimsunTimeRangeT = fcn_findValidTimeRangeUnixAndAimsun(names);
% 
% Author: Satya Prasad
% Created: 2020-05-28
% Updated: 2020-06-23
% 
% Revision history:
% 2020/05/29
%   - Added error check for inputs
% 2020/06/03
%   - Modified the output from global time range to unix time and bounds on 
%   system entrance time
% 2020/06/23
%   - Number of inputs are increased to two to have two types of queries.
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function unixTimeAndAimsunTimeRange = fcn_findValidTimeRangeUnixAndAimsun(names,name_road)
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
if ~(nargin == 1 || nargin == 2)
    error('fcn_findValidTimeRangeUnixAndAimsun: Incorrect number of input arguments');
end

%% Query for valid Unix Time and Bounds on System Entrance Time
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

if nargin == 1
    % SQL query to grab information of all trips
    sqlQuery_trip = ['SELECT id, date, description'...
                    ' FROM ' names.tableTrips];
else
    % SQL query to grab information of all trips on 'name_road'
    sqlQuery_trip = ['SELECT id, date, description'...
                    ' FROM ' names.tableTrips...
                    ' WHERE id IN (SELECT DISTINCT trips_id'...
                                 ' FROM ' names.tableTraffic...
                                 ' WHERE section_id IN (SELECT DISTINCT road_segment_id'...
                                                      ' FROM ' names.tableRoadSegment...
                                                      ' WHERE road_id = (SELECT id'...
                                                                       ' FROM ' names.tableRoad...
                                                                       ' WHERE name = ''' name_road ''')))'];
end
% grab trip information
tripInfo = fetch(DB.db_connection, sqlQuery_trip);

if ~isempty(tripInfo)
    % select date from the tripInfo
    timeDate = table2array(tripInfo(:,'date'));
    timeDate = datetime(timeDate);      % convert string to datetime object
    timeDate.TimeZone = 'America/New_York'; % set the time zone
    timeUnix = posixtime(timeDate);     % convert date to Unix Time
    % create a table containing Unix Time
    tableUnixTime = array2table(timeUnix,'VariableNames',{'UnixTime'});
    
    % initialize the variable to store bounds on Aimsun time
    timeAimsunBounds = ones(height(tableUnixTime),2);
    
    tripIds = table2array(tripInfo(:,'id'));    % trip IDs as an array
    for i = 1:length(tripIds)
        if nargin == 1
            % SQL query to grab min and max system entrance time for all trips
            sqlQuery_timeAimsun = ['SELECT MIN(system_entrance_time), MAX(system_entrance_time)'...
                                  ' FROM ' names.tableTraffic...
                                  ' WHERE trips_id = ' num2str(tripIds(i))];
        else
            % SQL query to grab min and max system entrance time for all 
            % trips on 'name_road'
            sqlQuery_timeAimsun = ['SELECT MIN(system_entrance_time), MAX(system_entrance_time)'...
                                  ' FROM ' names.tableTraffic...
                                  ' WHERE trips_id = ' num2str(tripIds(i))...
                                        ' AND section_id IN (SELECT DISTINCT road_segment_id'...
                                                           ' FROM ' names.tableRoadSegment...
                                                           ' WHERE road_id = (SELECT id'...
                                                                            ' FROM ' names.tableRoad...
                                                                            ' WHERE name = ''' name_road '''))'];
        end
        % grab min and max system entrance time
        result = fetch(DB.db_connection, sqlQuery_timeAimsun);
        % convert table to array
        timeAimsunBounds(i,:) = table2array(result);
    end
    % create a table containing bounds on Aimsun Time
    tableAimsunBounds = array2table(round(timeAimsunBounds,2),...
                        'VariableNames',{'aimsunLB','aimsunUB'});
    
    unixTimeAndAimsunTimeRange = [tripInfo(:,{'id','description'}), ...
                                  tableUnixTime, tableAimsunBounds];
else
    error('fcn_findValidTimeRangeUnixAndAimsun: No valid trip information');
end

% Disconnect from the database
DB.disconnect();

end