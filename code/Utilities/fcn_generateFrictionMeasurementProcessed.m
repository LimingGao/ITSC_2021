%%%%%%%%%%%%%%%%%%%%%  Function fcn_generateFrictionMeasurementProcessed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      given trip name and ENULLAref, then query raw friction measurement data from "friction_measurement" table
%      ,process it, and then insert them to "friction_measurement_processed" table 
%
% Input Variables:
%      names_trip      trip name, format:string
%      refLLA          ENU and LLA convert reference,format:3*1 array .
%
% Returned Results:
%
%
% Example:
% fcn_generateFrictionMeasurementProcessed('aimsun simluation 2020-05-13',[[-48.8766670000000,-123.393333000000,0]])
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
% Created Date:       2020-09-07
% Revisions:

% To do list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn_generateFrictionMeasurementProcessed(name_database,names_trip,refLLA)

%CONNECT TO  DATABASE
DB = Database(name_database);
% find the row numbers given a trip name
trip_id = fetch(DB.db_connection,['select id from trips where name = ''' names_trip '''']);
trip_id = trip_id{1,1};
rows_idRange = fetch(DB.db_connection,['select count(id),min(id),max(id) from friction_measurement where trips_id = ' num2str(trip_id)]);
row_nums = rows_idRange{1,1};

% process the data by a batch of 100,000 rows
insert_rows = 100000; % inssert 100,000 rows each loop
Split = ceil(row_nums/insert_rows); % inssert 10,000 rows each loop
tic
for i= 1:Split
    id_start = insert_rows*(i - 1) + rows_idRange{1,2};
    id_end = min(insert_rows*i - 1 + rows_idRange{1,2},rows_idRange{1,3});
    sql_friction_measurement = ['SELECT * FROM friction_measurement WHERE trips_id = ' num2str(trip_id) ' AND id >= ' num2str(id_start) ' AND id <= ' num2str(id_end)];
    friction_measurement_table = fetch(DB.db_connection,sql_friction_measurement);
    
    % {'id','trips_id','vehicle_id','vehicle_types_id','road_id','road_segment_id','road_lane_id','true_station_section','true_station_total','sensors_id','contact_point_east',
    % 'contact_point_north','contact_point_up','contact_point_latitude','contact_point_longitude','contact_point_altitude','contact_point_geography','friction_true','friction_measurement',
    % 'friction_measurement_noisy','noise_parameters','vehicle_model','tire_model','friction_model','unix_time_seconds','unix_time_nanoseconds','unix_time','timestamp','system_entrance_time',
    % 'aimsun_time','simulink_time','g_time','gps_week','gps_second','date_added'}
    
    %% front left tire
    % step1: parse contact point
    fl_contact_point = friction_measurement_table.fl_contact_point;
    fl_contact_point_cell = cell(length(fl_contact_point),3);
    for i= 1:length(fl_contact_point) %#ok<FXSET>
        ss = fl_contact_point{i};
        fl_contact_point_cell(i,:) = cell(ss.getArray())';
    end
    fl_contact_point_num = cell2mat(fl_contact_point_cell);
    % step2: to table
    friction_measurement_fl = friction_measurement_table(:,{'trips_id','vehicle_id','road_id','road_segment_id','unix_time_seconds','unix_time_nanoseconds','unix_time','timestamp', 'g_time'});
    
    friction_measurement_fl.contact_point_east = fl_contact_point_num(:,1);
    friction_measurement_fl.contact_point_north = fl_contact_point_num(:,2);
    friction_measurement_fl.contact_point_up = fl_contact_point_num(:,3);
    
    friction_measurement_fl.friction_true = friction_measurement_table.fl_friction_true;
    friction_measurement_fl.friction_measurement =friction_measurement_table.fl_friction_measurement;
    friction_measurement_fl.friction_measurement_noisy = friction_measurement_table.fl_friction_noisy;
    
    %%front right tire
    fr_contact_point = friction_measurement_table.fr_contact_point;
    fr_contact_point_cell = cell(length(fr_contact_point),3);
    for i= 1:length(fr_contact_point) %#ok<FXSET>
        ss = fr_contact_point{i};
        fr_contact_point_cell(i,:) = cell(ss.getArray())';
    end
    fr_contact_point_num = cell2mat(fr_contact_point_cell);
    % step2: to table
    friction_measurement_fr = friction_measurement_table(:,{'trips_id','vehicle_id','road_id','road_segment_id','unix_time_seconds','unix_time_nanoseconds','unix_time','timestamp', 'g_time'});
    
    friction_measurement_fr.contact_point_east = fr_contact_point_num(:,1);
    friction_measurement_fr.contact_point_north = fr_contact_point_num(:,2);
    friction_measurement_fr.contact_point_up = fr_contact_point_num(:,3);
    friction_measurement_fr.friction_true = friction_measurement_table.fr_friction_true;
    friction_measurement_fr.friction_measurement =friction_measurement_table.fr_friction_measurement;
    friction_measurement_fr.friction_measurement_noisy = friction_measurement_table.fr_friction_noisy;
    
    %%rear left tire
    % step1: parse contact point
    rl_contact_point = friction_measurement_table.rl_contact_point;
    rl_contact_point_cell = cell(length(rl_contact_point),3);
    for i= 1:length(rl_contact_point) %#ok<FXSET>
        ss = rl_contact_point{i};
        rl_contact_point_cell(i,:) = cell(ss.getArray())';
    end
    rl_contact_point_num = cell2mat(rl_contact_point_cell);
    % step2: to table
    friction_measurement_rl = friction_measurement_table(:,{'trips_id','vehicle_id','road_id','road_segment_id','unix_time_seconds','unix_time_nanoseconds','unix_time','timestamp', 'g_time'});
    
    friction_measurement_rl.contact_point_east = rl_contact_point_num(:,1);
    friction_measurement_rl.contact_point_north = rl_contact_point_num(:,2);
    friction_measurement_rl.contact_point_up = rl_contact_point_num(:,3);
    
    friction_measurement_rl.friction_true = friction_measurement_table.rl_friction_true;
    friction_measurement_rl.friction_measurement =friction_measurement_table.rl_friction_measurement;
    friction_measurement_rl.friction_measurement_noisy = friction_measurement_table.rl_friction_noisy;
    
    %%rear right tire
    % step1: parse contact point
    rr_contact_point = friction_measurement_table.rr_contact_point;
    rr_contact_point_cell = cell(length(rr_contact_point),3);
    for i= 1:length(rr_contact_point) %#ok<FXSET>
        ss = rr_contact_point{i};
        rr_contact_point_cell(i,:) = cell(ss.getArray())';
    end
    rr_contact_point_num = cell2mat(rr_contact_point_cell);
    % step2: to table
    friction_measurement_rr = friction_measurement_table(:,{'trips_id','vehicle_id','road_id','road_segment_id','unix_time_seconds','unix_time_nanoseconds','unix_time','timestamp', 'g_time'});
    
    friction_measurement_rr.contact_point_east = rr_contact_point_num(:,1);
    friction_measurement_rr.contact_point_north = rr_contact_point_num(:,2);
    friction_measurement_rr.contact_point_up = rr_contact_point_num(:,3);
    
    friction_measurement_rr.friction_true = friction_measurement_table.rr_friction_true;
    friction_measurement_rr.friction_measurement =friction_measurement_table.rr_friction_measurement;
    friction_measurement_rr.friction_measurement_noisy = friction_measurement_table.rr_friction_noisy;
    
    %% Concatenate them
    friction_measurement_processed =vertcat(friction_measurement_fl, friction_measurement_fr, friction_measurement_rl,friction_measurement_rr);
    % enu2lla
    [lat, lon, alt] = enu2geodetic(friction_measurement_processed.contact_point_east, friction_measurement_processed.contact_point_north, friction_measurement_processed.contact_point_up, refLLA(1),refLLA(2),refLLA(3),referenceEllipsoid('wgs84'));
    friction_measurement_processed.contact_point_latitude = lat;
    friction_measurement_processed.contact_point_longitude =lon ;
    friction_measurement_processed.contact_point_altitude = alt;
    
    friction_measurement_processed = movevars(friction_measurement_processed,{'unix_time_seconds','unix_time_nanoseconds','unix_time','timestamp', 'g_time'},'After','friction_measurement_noisy');
    
    %     friction_measurement_processed_insert.
    %     = friction_measurement_table(row_index_start:row_index_end,:);
    fprintf(1,' Insert: %.2f , %.2f. \n',id_start,id_end);
    if 1==1
        sqlwrite(DB.db_connection,'friction_measurement_processed',friction_measurement_processed)
    end
    toc
end
toc

%DISCONNECT TO DATABASE
DB.disconnect();


end
