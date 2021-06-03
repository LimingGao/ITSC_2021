%%%%%%%%%%%%%%%%%%%%%  Function fcn_frictionGridQueryBySection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      query friction data from database given road section id 
%
% Input Variables:
%      sectionID       road section id, format: array
%      path            the directory to save the query result(.mat file),
%                      this input is optional. 
%      method          (OPTIONAL) string flag - denotes the method used to query the data .
%
%                       method may be any of 'cluster', 'all'
%        
%                       method == 'cluster' --> query friction data from the cluster_friction_grid table.
%
%                       method == 'all' --> query friction data from the road_lane_grid table.
%
%                       DEFAULT: 'cluster'
% 
% Returned Results:
%      frictionGrid_matrix      the result in matrix format 
%      frictionGrid_table       the result in table format
%
% Example:
% [frictionGrid_matrix,frictionGrid_table] = fcn_frictionGrid(814,'.\DataFiles\')
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
% Created Date:       2020-05-28
% Revisions:
%           2020-05-29: 
%
% To do list: nargin
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frictionGrid_matrix,frictionGrid_table] = fcn_frictionGridQueryBySection(sectionID,path,varargin)

%step0: check input
% switch nargin
%     case 3
%         saveTofile = 1; % logical number to indicate if save the result to .mat file
%     case 2 
%         saveTofile = 1; % logical number to indicate if save the result to .mat file
%     case 1
%         saveTofile = 0;
%     otherwise
%         error('Incorrect number of input arguments');
% end

if nargin < 2
    error('Incorrect number of input arguments');
end

if strcmp(path,'None') || strcmp(path,'none')
    saveTofile  = 0; % flag
else
    saveTofile  = 1; % flag
    if numel(path) > 0 && path(end) ~= '/'
        path = [path '/'];
    end
end

method = 'cluster'; % nargin==2

% are there any other arguments?
if nargin > 2
  % there are. check the last argument. Is it a string?
  if ischar(varargin{end})
    method = varargin{end};
%     varargin(end) = [];
  else 
      method = 'cluster';
  end

end

if method(1) == 'c'
    table_name = 'cluster_friction_grid';
    fields = 'id,road_lane_id,road_lane_grid_id,position_x,position_y,position_z,latitude,longitude,altitude,friction_coefficient,station';

    fileName_suffix = '_cluster_';  
elseif method(1) == 'a'
    table_name = 'road_lane_grid';  
    fields = 'id,road_lane_id,lane_surface_condition_id,position_x,position_y,position_z,latitude,longitude,altitude,friction_coefficient';

    fileName_suffix = '_';
end

%% step1: CONNECT TO  DATABASE ----------------- %
% choose different da name to connect to them
database_name = 'nsf_roadtraffic_friction_v2';

% create a instance
DB = Database(database_name);

%% step2: query the friction grid from DB ----------------- %
sectionID = sectionID(:);
disp(['Querying frictionGrid data of road section ' num2str(sectionID') ' from ' database_name ' ...']);
tic

where = strjoin(cellstr(num2str(sectionID)), ', '); % char
where = cat(2,'road_segment_id in (', where, ')'); % 1*1 cell

sql_friction =[ 'select ' fields ' from ' table_name ' where road_lane_id in '...
                '(select road_lane_id from road_segment_lane where ' where ')'];
% 
% sql_friction =[ 'select ' fields ' from road_lane_grid where road_lane_id in '...
%                 '(select road_lane_id from road_segment_lane where road_segment_id in ( ' num2str(sectionID) '))'];
frictionGrid_table = fetch(DB.db_connection,sql_friction);

frictionGrid_matrix = table2array(frictionGrid_table);
toc
disp(['Query done! ' num2str(height(frictionGrid_table)) ' rows frictionGrid data quried from ' database_name ' !']);

%%
if saveTofile
    save(strcat(path,'frictionGrid_table_section',fileName_suffix,strjoin(cellstr(num2str(sectionID)),'_'),'.mat'),'frictionGrid_table')
    save(strcat(path,'frictionGrid_matrix_section',fileName_suffix,strjoin(cellstr(num2str(sectionID)),'_'),'.mat')','frictionGrid_matrix')
    
    disp('frictionGrid data has been saved to /DataFiles.' );
end

%% disconnect
DB.disconnect();

end


