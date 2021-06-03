%%%%%%%%%%%%%%%%%%%%%  Function fcn_laneTypesQueryByLaneIds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      query lane types from database given road lane Ids
%
% Input Variables:
%      laneID          road lane ids, format:array
% 
% Returned Results:
%      laneTypes       the lane type of each lane, format:table
%
% Example:
% [laneTypes] = fcn_laneTypesQueryByLaneIds(20101)
% 
% Processing Flow:
%
% Restrictions/Notes:
%
% The following functions are called:
%      none
%
% Author:             Liming Gao
% Created Date:       2020-05-28
% Revisions:
%           2020-05-29: 
%
% To do list: nargin
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [laneTypes] = fcn_laneTypesQueryByLaneIds(laneID)

%step0: check input

if nargin < 1
    error('Incorrect number of input arguments');
end

%% step1: CONNECT TO DATABASE ----------------- %
% choose different da name to connect to them
database_name = 'nsf_roadtraffic_friction_v2';

% create a instance
DB = Database(database_name);

DB.db_connection.AutoCommit
%% step2: query the lane properties from DB ----------------- %
laneID = laneID(:);

disp(['Querying lane types data of road lanes ' num2str(laneID') ' from ' database_name ' ...']);
tic

table_name = 'road_lane_types as R, lane_types as L'; % target table 
fields = 'road_lane_id,lane_types_id,name,access_vehiceltype,rule,description';
where = strjoin(cellstr(num2str(laneID)), ', '); % char
where = cat(2,'road_lane_id in (', where, ')'); % 1*1 cell

sql_laneTypes =[ 'select ' fields ...
                 ' from ' table_name ...
                 ' where R.lane_types_id = L.id AND ' where ];
    
laneTypes = fetch(DB.db_connection,sql_laneTypes);

toc
disp(['Query done! ' num2str(height(laneTypes)) ' rows data quried from ' database_name ' !']);

%% disconnect
DB.disconnect();

end


