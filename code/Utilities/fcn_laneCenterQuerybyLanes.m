%%%%%%%%%%%%%%%%%%%%%  Function fcn_laneCenterQuerybyLanes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      query lane center node from database given road section id 
%
% Input Variables:
%      laneID          road lanes id, format:array
%      path            the directory to save the query result(.mat file),
%                      this input is optional. 
%      
% Returned Results:
%      lanesCenter_matrix      the result in matrix format 
%      lanesCenter_table       the result in table format
%      
% Example:
% [lanesCenter_matrix,lanesCenter_table] = fcn_laneCenterQuerybyLanes(20101,'')
% 
% Processing Flow:
%
%
% Restrictions/Notes:
%
% The following functions are called:
%      none
%
% Author:             Liming Gao
% Created Date:       2020-06-03
% Revisions:
%           2020-06-05: 
%
% To do list:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lanesCenter_matrix,lanesCenter_table] = fcn_laneCenterQuerybyLanes(laneID,path)

%%step0: check input
switch nargin
    case 2
        if strcmp(path,'None') || strcmp(path,'none')
            saveTofile  = 0; % flag
        else
            saveTofile  = 1; % flag
            if numel(path) > 0 && path(end) ~= '/'
                path = [path '/'];
            end
        end
    case 1
        saveTofile = 0;
    otherwise
        error('Incorrect number of input arguments');
end

% Notes: for future 
% if (nargin==1)
%    z = sqrt(x);
% elseif (nargin==2)
%    z = sqrt(x+y);
% end


%% step1: CONNECT TO  DATABASE ----------------- %
% choose different da name to connect to them
database_name = 'nsf_roadtraffic_friction_v2';

% create a instance
DB = Database(database_name);

%% step2: query the lane Center node from DB ----------------- %
laneID = laneID(:);
disp(['Querying lane Center node data of road lanes ' num2str(laneID') ' from ' database_name ' ...']);
tic
%query lane node (lane centerline)
fields_lane_node =  {'road_lane_id','lane_node_id','center_x','center_y','center_z','latitude','longitude','altitude','lane_width','lane_border_left','lane_border_right','bank','grade','curvature','yaw','station'};
fields_lane_node = strjoin(fields_lane_node,', ');

where = strjoin(cellstr(num2str(laneID)), ', '); % char
where = cat(2,'road_lane_id in (', where, ')'); % 1*1 cell

sql_lane_node =[ 'select ' fields_lane_node ' from lane_node where ' where];
lanesCenter_table = fetch(DB.db_connection,sql_lane_node);

lanesCenter_matrix = table2array(lanesCenter_table);

toc
disp(['Query done! ' num2str(height(lanesCenter_table)) ' rows lane node data quried from ' database_name ' !']);

%%
if saveTofile
    save(strcat(path,'lanesCenter_table_lanes',strjoin(cellstr(num2str(laneID)),'_'),'.mat'),'lanesCenter_table')
    save(strcat(path,'lanesCenter_matrix_lanes',strjoin(cellstr(num2str(laneID)),'_'),'.mat')','lanesCenter_matrix')
    
    disp('lanesCenter data has been saved to /DataFiles.' );
end

%% disconnect
DB.disconnect();

end
