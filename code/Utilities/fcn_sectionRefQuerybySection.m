%%%%%%%%%%%%%%%%%%%%%  Function fcn_sectionRefQuerybySection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      query road section reference curve given road section id 
%
% Input Variables:
%      sectionID       road section id, format:array
%      path            the directory to save the query result(.mat file),
%                      this input is optional. 
%      
% Returned Results:
%      sectionRef_matrix      the result in matrix format 
%      sectionRef_table       the result in table format
%      roadLaneID              the lanes id in the sections 
% Example:
% [lanesCenter_matrix,lanesCenter_table] = fcn_laneCenterQuerybySection(814,'')
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

function [sectionRef_matrix,sectionRef_table] = fcn_sectionRefQuerybySection(sectionID,path)

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
sectionID = sectionID(:);
disp(['Querying road section reference node data of road section ' num2str(sectionID') ' from ' database_name ' ...']);
tic
%query lane node (lane centerline)
fields_section_ref =  {'road_segment_id','road_segment_reference_id','east','north','up','latitude','longitude','altitude','bank','grade','curvature','yaw','station'};
fields_section_ref = strjoin(fields_section_ref,', ');

where = strjoin(cellstr(num2str(sectionID)), ', '); % char
where = cat(2,'road_segment_id in (', where, ')'); % 1*1 cell

sql_section_ref =[ 'select ' fields_section_ref ' from road_segment_reference where ' where];
sectionRef_table = fetch(DB.db_connection,sql_section_ref);

sectionRef_matrix = table2array(sectionRef_table);

toc
disp(['Query done! ' num2str(height(sectionRef_table)) ' rows segemnt reference data quried from ' database_name ' !']);

%%
if saveTofile
    save(strcat(path,'sectionRef_table_section',strjoin(cellstr(num2str(sectionID)),'_'),'.mat'),'sectionRef_table')
    save(strcat(path,'sectionRef_matrix_section',strjoin(cellstr(num2str(sectionID)),'_'),'.mat')','sectionRef_matrix')
    
    disp('section Reference data has been saved to /DataFiles.' );
end

%% disconnect
DB.disconnect();

end
