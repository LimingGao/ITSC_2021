%%%%%%%%%%%%%%%%%%%%%  Function fcn_ENUrefPointQuery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      query ENU reference point given road name  
%
% Input Variables:
%      RoadName        road name, string
%      
% Returned Results:
%      referencePoint           ENU reference point,1*3 array [L L A]
% 
% Example:
% RefPoint = fcn_ENUrefPointQuery('Onekm_straight_highway')
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
% To do list:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RefPoint = fcn_ENUrefPointQuery(RoadName)

%%step0: check input, Note:more edit for the the input of road id 
switch nargin
    case 2
        
    case 1
        flag_queryByRoadName = 1;
    otherwise
        error('Incorrect number of input arguments');
end


%% step1: CONNECT TO  DATABASE ----------------- %
% choose different da name to connect to them
database_name = 'nsf_roadtraffic_friction_v2';

% create a instance
DB = Database(database_name);

%%show tables
tables = DB.ShowTables(); %#ok<NASGU>

%% step2: query the referecne lla from DB ----------------- %
if isequal(flag_queryByRoadName,1)
    fprintf(1,'\n')
    disp(['Querying ENU reference point for road: ' RoadName]);
    
    fields = 'id,name,latitude,longitude,altitude';
    
    sql_enuRef =[ 'select ' fields ' from enu_reference where id in '...
        '(select enu_reference_id from road where name = ' char(39) RoadName char(39) ')'];
    enuRef_table = fetch(DB.db_connection,sql_enuRef);
    
    disp(['Query done! The road ' RoadName ' use ' enuRef_table.name{1} ' ENU reference Point!']);
    fprintf(1,'\n')
end
%%result 
RefPoint = [enuRef_table.latitude, enuRef_table.longitude, enuRef_table.altitude];

%% disconnect
DB.disconnect();

end


