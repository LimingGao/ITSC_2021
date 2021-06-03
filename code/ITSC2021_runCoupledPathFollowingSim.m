%%%%%%%%%%%%  Script ITSC2021_runCoupledPathFollowingSim.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used for the ITSC2021 conference
% Function Purpose: Run Simulink ITSC2021_CoupledPathFollowingSim_discrete.slx
% Matlab work Path: ~\GitHub\forgetfulDBs
% Author:       Liming, Satya, Dr. Beal, Dr.Brennan
% Created Date: 2020-05-15
%
% Revisions:
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: Prepare the workspace
% Clear workspace, figures, and console
clear all;  %#ok<CLALL> % Clears workspace to remove any old variables
close all;  % Close all open figures
clc;        % Clear console space
%%
addpath('./Utilities/');    % all the functions and wrapper class
addpath('./DataFiles/');    % all the .mat data files
addpath('./Generate_longitudinal_velocity_profile/');    % all the .mat data files
addpath(genpath('./Figures/'));    % all the image files 

dir.datafiles = ['.' filesep 'DataFiles' filesep]; % all the .mat data files

% load dirsed path
load('path_desired.mat','path'); 
% remove the Duplicate data

[~,unique_station_indices] = unique(path.station,'stable');
subfieldNames = fieldnames(path); % Grab all the subfields
for i_subField = 1:length(subfieldNames)
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};
        path.(subFieldName) = path.(subFieldName)(unique_station_indices); % repalce with unique data    
end
          
% 'names' is a structure containing name(s) of database, tables/relations
% ---------------------------------------------------

% ----------------------------------------------------
% name of the database
names.database = 'nsf_roadtraffic_friction_v2';
% name of the table that contains the relation between road name and road id
names.tableRoad = 'road';
% name of the table that contains the relation between road id and
% section/segment ids
names.tableRoadSegment = 'road_road_segment';
% name of the table containing traffic data
names.tableTraffic = 'road_traffic_procsessed';

%referenceEllipsoid
spheroid = referenceEllipsoid('wgs84');

% flag triggers
flag.dbQuery = false; % set to 'true' to query new data from DB. Otherwise Onekm_curved_highway default file will be loaded.
flag.simByvehicleID = false; % only simulate the vehicles with specific ids
flag.dbInsert = false; % set to 'true' to insert friction measurement to database 
flag.saveTrajectory = false; % set to 'true' to save vehicle trajectory as a mat file
flag.saveFricMeasure = false; % set to 'true' to save friction_measurement as a mat file
flag.verbose = true; % flag to print to the console
flag.doDebug = true; % set to 'true' to print trajectory information to command window
flag.plot = true; % flag for plot
flag.noise = true; % flag to add noise to sensor measurements


% swerve parameters for adding two sine waves
flag.swerve        = true;  % 'true' to add distortion and 'false' to not
flag.lateralOffset = true;  % 'true' to add random offset of the trajectory and 'false' to not
% amplitude(s) of sine distortion in [meters]
% perpendicular to center of the aimsun vehicle trajectory
swerve.amplitude  = [0.1, 0.4]; %[meters] [0.3,0.1]
% time period(s) of sine distortion in [seconds]
swerve.timePeriod = [37, 4]; %[seconds][60,5]
% maximum fixed lateral offset to a vehicle trajectory
swerve.MaxLateralOffset = 0.1; %[meters]

% time series parameters
flag.timeSeries = false;% 'true' to generate time series and 'false' to not
timeLookAhead   = 0;    % look ahead time in [seconds]

% controller status
flag.controller = false; % set to 'true' to active longitudinal control
flag.speedUp = false; % set to 'true' to speedup the vehicle velocity 
velocity_factor = 1; % the scale to speed up the ego vehicle
flag.laneCenterTraj = true; % set to 'true' if you want the vehicle follow exactly one lane center line
flag.speed_preview =false;


% field names of vehicle trajectory before being used in the simulation
fieldsPreSimulationVehicleTrajectory.vehId      = 1;
fieldsPreSimulationVehicleTrajectory.vehType    = 2;
fieldsPreSimulationVehicleTrajectory.vehLength  = 3;
fieldsPreSimulationVehicleTrajectory.vehWidth   = 4;
fieldsPreSimulationVehicleTrajectory.sectiondId = 5;
fieldsPreSimulationVehicleTrajectory.junctionId = 6;
fieldsPreSimulationVehicleTrajectory.laneNum    = 7;
fieldsPreSimulationVehicleTrajectory.direction  = 8;
fieldsPreSimulationVehicleTrajectory.lat        = 9;
fieldsPreSimulationVehicleTrajectory.lon        = 10;
fieldsPreSimulationVehicleTrajectory.alt        = 11;
fieldsPreSimulationVehicleTrajectory.east       = 12;
fieldsPreSimulationVehicleTrajectory.north      = 13;
fieldsPreSimulationVehicleTrajectory.up         = 14;
fieldsPreSimulationVehicleTrajectory.yaw        = 15;
fieldsPreSimulationVehicleTrajectory.speed      = 16;
fieldsPreSimulationVehicleTrajectory.station    = 18;
fieldsPreSimulationVehicleTrajectory.aimsunTime = 19;
fieldsPreSimulationVehicleTrajectory.globalTime = 20;

% reference LLA for ENU to LLA transformation
if flag.dbQuery
    refLLA = fcn_ENUrefPointQuery(names.road); % query corresponding reference point from dababase
else
    % NEMO Point, Centre of pacific ocean, is used as reference for virtual
    % simulations
    refLLA = [-48.876667, -123.393333, 0];
end

%% step 2: Get/Check the validity of Unix Time and Bounds on System Entrance Time
% Unix Time is the number of seconds that have elapsed since 1970-01-01 00:00:00.
% It is computed from datetime using the function 'posixtime'. It is in [seconds]

% System Entrance Time is the time at which a vehicle enters the network.
% It is not wall time but the time measured by AIMSUN or time elapsed in
% the simulation. It is in [seconds]

% Global Time is Unix Time plus System Entrance Time

% unixTimeAndAimsunTimeRange: It's a Nx5 table. First and second columns
% contian trip id and trip description. Third column contains Unix Time.
% Fourth and fifth columns contains lower and upper bounds on System
% Entrance Time respectively. Output contains information of all trips if
% number of inputs is one. Output contains information of all trips on
% 'names.road' if number of inputs is two.
if flag.dbQuery
    unixTimeAndAimsunTimeRange = fcn_findValidTimeRangeUnixAndAimsun(names,names.road);

    indexTripOfInterest = 1;
    % UNIX Time at which simulation is started, [seconds]
    time.Unix = unixTimeAndAimsunTimeRange(indexTripOfInterest,'UnixTime');
    % lower bound on the System Entrance Time, [seconds]
    time.AimsunLB = unixTimeAndAimsunTimeRange(indexTripOfInterest,'aimsunLB');
    % upper bound on the System Entrance Time, [seconds]
    time.AimsunUB = unixTimeAndAimsunTimeRange(indexTripOfInterest,'aimsunUB');
else
    indexTripOfInterest = 0;
    time.Unix = 0;
    time.AimsunLB = 0;
    time.AimsunUB = 0; 
end
    
% if flag.verbose
%     % Print 'unixTimeAndAimsunTimeRange' to command window
%     fprintf(1,'\nUnix Time(s)   Lower Bound(s)   Upper Bound(s) \n');
%     fprintf(1,[repmat('%d       %.2f             %.2f \n', 1, ...
%         size(unixTimeAndAimsunTimeRange,2)) '\n'], unixTimeAndAimsunTimeRange');
% end % NOTE: Ends flag.verbose
%
% if flag.doDebug
%     time.Unix = 1589072657; % UNIX Time at which simulation is started, [seconds]
%     time.AimsunLB = 6.61;   % lower bound on the System Entrance Time, [seconds]
%     time.AimsunUB = 1795.61;% upper bound on the System Entrance Time, [seconds]
%
%     % checks the validity of Unix Time, lower and upper bounds on System
%     % Entrance Time
%     fcn_checkValidityUnixAndAimsunRange(unixTimeAndAimsunTimeRange,time)
% end % NOTE: Ends flag.doDebug

%% Step 3: Get all the valid Vehicle ID, Global Time combinations lying
% within specified bounds of global time defined by Unix Time and System
% Entrance Time
% vehiclePopulation_ID_GlobalTime: Nx2 vector, consisting of unique
% (vehicle_id, global_time) pairs. N is number of unique trajectories.

% SectionId_VehID_GlobalTime: Mx3 vector, consisting of unique
% (section_id, vehicle_id, global_time) pairs. M is number of unique
% encounters of all trajectories with sections in a road.

% listOfSections: Lx1 vector,consisting of unique section_ids

% just query one trajectory if flag.controller is true
if ~flag.controller && flag.dbQuery
    [~,SectionId_VehID_GlobalTime,listOfSections] = ...
        fcn_findValidVehicleIdandGlobalTime(time.Unix,time.AimsunLB,...
        time.AimsunUB,names,names.road);
elseif flag.controller && flag.dbQuery
    %SectionId_VehID_GlobalTime = [201,23,1589416895.82000];
    SectionId_VehID_GlobalTime = [201,51,1589416967.21000];
    listOfSections = 201;
    
elseif ~flag.dbQuery % load the default trajectory
    SectionId_VehID_GlobalTime = [201 593 1589418271]; %load default files 
    listOfSections = 201;
end

% hardcode for dubugging to a single vehicle trajectory
% listOfSections = 201;   % 
% indexTrajectoryOfInterest = 24; row number in 'vehiclePopulation_ID_GlobalTime'

    %% step 4: query the friction grid and lane center data or load from .mat files
    % friction_section_id: sections of a road that are included in the vehicle trajectory
    % frictionGrid_table: [grid location in ENU, friction coefficient] in table format
    % loading friction data from the database depending of the sections part of the vehicle trajectory
    % find road section from SectionId_VehID_GlobalTime given a vehicle_id
    if flag.dbQuery
        % query friction grid data given road sections
        if flag.saveTrajectory
            data_path = dir.datafiles;
        else
            data_path = 'None';
        end
        % query road gird and friction values
        [~,frictionGrid_table] = fcn_frictionGridQueryBySection(listOfSections(indexSectionOfInterest),data_path,'all');
        [~,frictionGrid_Cluster_table] = fcn_frictionGridQueryBySection(listOfSections(indexSectionOfInterest),data_path,'cluster');
        
        % query lane center data given road sections
        [~,lanesCenter_table] = fcn_laneCenterQuerybySection(listOfSections(indexSectionOfInterest),data_path);
        
        % query road segment reference data 
        [~,sectionRef_table] = fcn_sectionRefQuerybySection(listOfSections(indexSectionOfInterest),data_path);
        
        % check queried data
        if flag.doDebug
            fcn_plotLoadFrictionGrid(frictionGrid_table,lanesCenter_table,frictionGrid_Cluster_table)
        end
    else
        % load default data
        % load('frictionGrid_table_section814.mat','frictionGrid_table') % load friction grid data
        load('frictionGrid_table_section_cluster_201.mat','frictionGrid_table')  % load friction grid cluster data
        frictionGrid_Cluster_table = frictionGrid_table;
        
        
        frictionGrid_Cluster_table = frictionGrid_Cluster_table(1:length(path.x),:);
        frictionGrid_Cluster_table.id = (1:1:length(path.x))';
        frictionGrid_Cluster_table.position_x = path.x;
        frictionGrid_Cluster_table.position_y = path.y;
        frictionGrid_Cluster_table.position_z = zeros(size(path.y));
        
        frictionGrid_Cluster_table.friction_coefficient = path.friction_true;
        frictionGrid_Cluster_table.station = path.station;
        
        frictionGrid_Cluster_table.road_lane_grid_id = zeros(size(path.y));
        frictionGrid_Cluster_table.latitude = zeros(size(path.y));
        frictionGrid_Cluster_table.longitude = zeros(size(path.y));
        frictionGrid_Cluster_table.altitude = zeros(size(path.y));
        
        frictionGrid_table = frictionGrid_Cluster_table;
        
%         load('frictionGrid_table_section_201.mat','frictionGrid_table') % load friction grid data
        load('lanesCenter_table_section201.mat','lanesCenter_table') % load lane center data
        
        load('sectionRef_table_section201.mat','sectionRef_table') % load lane center data
        
        % check loaded data
        if flag.doDebug
%             fcn_plotLoadFrictionGrid(frictionGrid_table,lanesCenter_table,frictionGrid_Cluster_table)
        end
    end
    
    % prepare data for ST Friction DB Query block
    % simulink model need this friction_grid constant[xEast, yNorth, zUp, friction_coefficient]
    
    % friction_grid with cluster
    friction_grid_cluster = [frictionGrid_Cluster_table.position_x, frictionGrid_Cluster_table.position_y, frictionGrid_Cluster_table.position_z, frictionGrid_Cluster_table.friction_coefficient]; % cluster grid position and friction
    %index_change = find(friction_grid_cluster(:,4)==0.1); % Note:need to remove
    %friction_grid_cluster(index_change,4)=0.12;
    
    % prepare data for ST Road Geometry DB Query block 
    % road_geometry
    road_geometry_radius =1./sectionRef_table.curvature;% smoothdata(1./sectionRef_table.curvature,'gaussian',180);
    road_geometry = [sectionRef_table.east, sectionRef_table.north, sectionRef_table.up, sectionRef_table.grade, sectionRef_table.bank, road_geometry_radius, sectionRef_table.yaw]; % 
       
    % Note: this is just for plotting and debuging
    Vehicle_Geometry = [1.5 ;2.25; 2; 0.5]; % a,b,T,h, a/b = 0.4/0.6,simulink might use this,
    Wheel_radius = 0.32; %[meters] Wheel radius of vehicle
    
    % find vehiclePopulation_ID_GlobalTime by section
%     vehiclePopulation_ID_GlobalTime = ...
%         SectionId_VehID_GlobalTime((SectionId_VehID_GlobalTime(:,1)==...
%         listOfSections(indexSectionOfInterest)),2:3);
    
    
    if flag.simByvehicleID % simulate by vehicle id 
        %[row,col] = find(vehiclePopulation_ID_GlobalTime(:,1)<=120 & vehiclePopulation_ID_GlobalTime(:,1)>=100); % vehiclePopulation_ID_GlobalTime(:,1) indicates the vehicle ids
        %[row,col] = find(vehiclePopulation_ID_GlobalTime(:,1)<=100 & vehiclePopulation_ID_GlobalTime(:,1)>=1); % simulated id from 1 to 100 of 'Tilted_T' road on 2021-01-23 
        [row,col] = find(vehiclePopulation_ID_GlobalTime(:,1)<=200 & vehiclePopulation_ID_GlobalTime(:,1)>=101); % simulated id from 101 to 200 of 'Tilted_T' road on 2021-01-24
        %[row,col] = find(vehiclePopulation_ID_GlobalTime(:,1)==100); % vehiclePopulation_ID_GlobalTime(:,1) indicates the vehicle ids
        
        TrajectoryOfInterest = row(:)'; % index of specific vehicle ids 
    else % simulate all vehicles
        TrajectoryOfInterest = [1]; % index of all vehicle ids 
    end %// end of flag.simByvehicleID
    
    for indexTrajectoryOfInterest = TrajectoryOfInterest % inner loop of trajectory,eg. 12
        %% step 5: load a vehicle's trajectory from database or a mat file
        % queryVehicleTrajectory: Nx20 matrix, ontaining all the attributes
        % defined by trajectory_attributes(within the function) sorted in
        % the order of aimsun_time
        if flag.dbQuery
            % query vehicle trajectory from the database
            queryVehicleTrajectory = fcn_queryVehicleTrajectory(...
                vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,1), ...
                vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,2), ...
                names, listOfSections(indexSectionOfInterest));
            % save the trajectory as a mat file, Just in case to dubug it offline
            if flag.saveTrajectory
                save([dir.datafiles, 'vehicle_trajectory_', ...
                    num2str(vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,1)), '_', ...
                    num2str(round(vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,2))), ...
                    '.mat'], 'queryVehicleTrajectory');
            end %  Ends flag.saveTrajectory
        else
            % load queryVehicleTrajectory data, Just in case to dubug it offline
            load('vehicle_trajectory_26_1589072708.mat','queryVehicleTrajectory'); 
        end % Ends flag.dbQuery
        
        % check loaded trajectory
        if flag.doDebug
%              fcn_plotTrajectory(queryVehicleTrajectory); 
        end
         
        if flag.verbose
            % print the duration of vehicle trajectory
            fprintf(1,'\nDuration of vehicle - %d trajectory is %.2f seconds \n', ...
                queryVehicleTrajectory(1,1), ...
                max(queryVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.aimsunTime))-...
                min(queryVehicleTrajectory(:,fieldsPreSimulationVehicleTrajectory.aimsunTime)));
            % print the start and end point of vehicle trajectory
            fprintf(1,'Trajectory of vehicle - %d begins at (%.2f, %.2f, %.2f) and ends at (%.2f, %.2f, %.2f) \n', ...
                queryVehicleTrajectory(1,1), ...
                queryVehicleTrajectory(1,fieldsPreSimulationVehicleTrajectory.east), ...
                queryVehicleTrajectory(1,fieldsPreSimulationVehicleTrajectory.north), ...
                queryVehicleTrajectory(1,fieldsPreSimulationVehicleTrajectory.up), ...
                queryVehicleTrajectory(end,fieldsPreSimulationVehicleTrajectory.east), ...
                queryVehicleTrajectory(end,fieldsPreSimulationVehicleTrajectory.north), ...
                queryVehicleTrajectory(end,fieldsPreSimulationVehicleTrajectory.up));
        end % Ends flag.verbose
        
        %% step 6: Add swerving to the queryVehicleTrajectory
        % fcn_swervedVehicleTrajectory takes in queryVehicleTrajectory and adds two
        % sine waves to it based on swerve.amplitude, swerve.time_period
        % fixed lateral offset to the vehicle trajectory
        if flag.lateralOffset
           swerve.lateralOffset = swerve.MaxLateralOffset*(2*rand-1); % add random lateral offset
        else
            swerve.lateralOffset = 0;
        end
        swervedVehicleTrajectory = queryVehicleTrajectory;
%         = fcn_swervedVehicleTrajectory(queryVehicleTrajectory,...
%             swerve.amplitude, swerve.timePeriod, swerve.lateralOffset, refLLA, ...
%             fieldsPreSimulationVehicleTrajectory, flag.swerve); % lateral offset works even flag.swerve is false
        % Notes: in this function, yaw is changed to atan2(diff(new_north), diff(new_east));
        %% step 7: create time series data for vehicle trajectory
        % if flag.timeSeries is 'true' -> we are doing time based vehicle simulations
        % otherwise space based vehicle simulations
        if  flag.timeSeries
            % trajectory is time based
            [timeSeriesVehicleTrajectory, sim_duration_from_aimsun_time] = ...
                fcn_timeSeriesVehicleTrajectory(swervedVehicleTrajectory, ...
                fieldsPreSimulationVehicleTrajectory, timeLookAhead);
            % input to the simulink model
            inputTrajectory = timeSeriesVehicleTrajectory.Data;
        else
            % trajectory is space based input to the simulink model
            inputTrajectory = swervedVehicleTrajectory; 
            
            % due to the yaw jumps of trajectory ffrom aimsun, I use Lane center line as the reference trajectory

            if flag.laneCenterTraj

                inputTrajectory = repmat(swervedVehicleTrajectory,5,1);

                %position interplation
                nb_query = length(inputTrajectory);
                pt_en = fcn_interparc(nb_query,path.x,path.y,'pchip');
                pt_u = zeros(nb_query,1);
                % calcualte station
                inputTrajectory_Diff = [[0 0]; diff(pt_en(:,1:2))];
                inputTrajectory_Station = cumsum(sqrt(sum(inputTrajectory_Diff.^2,2)));
                % calcualte yaw
                pt_yaw_diff = atan2(diff(pt_en(:,2)), diff(pt_en(:,1)));
                pt_yaw = [pt_yaw_diff(1); pt_yaw_diff];
                
                % velocity interplation
                if flag.speed_preview
                    inputTrajectory_velocity = interp1(path.station,path.velocity_safe,inputTrajectory_Station,'linear','extrap');
                else
                    inputTrajectory_velocity = interp1(path.station,path.velocity_unsafe,inputTrajectory_Station,'linear','extrap');
              
                end
                    % update inputTrajectory
                inputTrajectory(:,[12:18]) = [pt_en pt_u pt_yaw inputTrajectory_velocity inputTrajectory_Station inputTrajectory_Station];
                % aimsun time 
                inputTrajectory(:,19) = (0.1:0.1:0.1*length(inputTrajectory))';
           
            end
            
            % change the speed
            if flag.speedUp
                ego_vehicle_speed = velocity_factor*inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.speed);
                inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.speed) = ego_vehicle_speed;
                
                % interpolate the trajectory to smooth because the trajectory from
                % aimsun is not smooth, which is connected by linear
                % segment.(But we ca not do this, this is helpful for the more fine preview location )
                %{
                [~,unique_station_indices] = unique(inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.station),'stable');
                inputTrajectory  = inputTrajectory(unique_station_indices,:);
                
                nb_query = length(inputTrajectory);
                pt_enu = fcn_interparc(nb_query,inputTrajectory(:,12),inputTrajectory(:,13),inputTrajectory(:,14),'pchip');
                
                pt_yaw_diff = atan2(diff(pt_enu(:,2)), diff(pt_enu(:,1)));
                pt_yaw = [pt_yaw_diff(1); pt_yaw_diff];
                pt_yaw_smooth = smoothdata(pt_yaw,'gaussian',38);
                
                inputTrajectory_Diff = [[0 0 0]; diff(pt_enu(:,1:3))];
                inputTrajectory_Station = cumsum(sqrt(sum(inputTrajectory_Diff.^2,2)));
                % need more edit
                %inputTrajectory(:,[12:15,18]) = [pt_enu pt_yaw inputTrajectory_Station]; % replace with interplation data
                
                
                figure(12126)
                plot(inputTrajectory_Station,pt_yaw,'b')
                hold on 
                plot(inputTrajectory_Station,smoothdata(pt_yaw,'gaussian',40),'r')
                %R = smoothdata(R,'gaussian',500); % smooth the result
                
                lanesCenter_yaw_diff = atan2(diff(lanesCenter_table.center_y), diff(lanesCenter_table.center_x));
                lanesCenter_yaw = [lanesCenter_yaw_diff(1); lanesCenter_yaw_diff];
                %figure(127126)
                plot(lanesCenter_table.station(1:9992),lanesCenter_yaw(1:9992),'g')
                grid on
                %}
                % duration of the simulation
                sim_duration_from_aimsun_time = round((max(inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.aimsunTime)) - ...
                    min(inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.aimsunTime)))*0.65);
            else
                % duration of the simulation
                sim_duration_from_aimsun_time = ...
                    max(inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.aimsunTime)) - ...
                    min(inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.aimsunTime));
            end
 
        end % NOTE: Ends flag.timeSeries
        
                % check loaded trajectory
        if flag.doDebug
%              fcn_plotTrajectory(queryVehicleTrajectory); 
        end
        
        %% step 8: run the simulation
        % set parameters in Friction query block
        mdlparam.frictionGridSize = size(friction_grid_cluster);  % size of friction grid
        mdlparam.refLLA = refLLA; % reference LLA for ENU coordinates
        
        % set parameters for road geometry block
        mdlparam.road_geometry = size(road_geometry);  % size of friction grid road_geometry
        mdlparam.road_geometry_preview_distance = 120; % distance of preview
        
        % set parameters in Trajectory query block   
        mdlparam.trajectorySize = size(inputTrajectory);    % size of input trajectory
        mdlparam.flagFullWndow  = false; % parameter that defines the search window for trajectory query
        % parameter to decide the size of the search window by using the maximum
        % distance that can be traversed by a vehicle moving with velocity 40m/s in 0.1s
        mdlparam.halfSearchDistance = 40*0.1; %[meters]
        
        mdlparam.lookAheadDistance  = 20; % look ahead distance[meters]
        mdlParam.horizonAVAR = 2048; % horizon size to calculate AVAR[samples]
        mdlParam.transformationWindow = 50; % window size to estimate friction transformation
        if flag.noise
            mdlParam.sampling_time = 0.01; % [seconds]
            mdlParam.sampling_time_accel = 0.01;
        else
            mdlParam.sampling_time = 0.0; % [seconds]
            mdlParam.sampling_time_accel = 0.0;
        end
        mdlParam.contactPatchLength = 0.15; % length of contact patch [meters]
        mdlParam.frictionCoefficientRatio = 1;
        % initialization
        initial.east  = inputTrajectory(1,fieldsPreSimulationVehicleTrajectory.east); %[meters]
        initial.north = inputTrajectory(1,fieldsPreSimulationVehicleTrajectory.north); %[meters]
        
        [~,initialAhead.index] = min(abs(inputTrajectory(:,fieldsPreSimulationVehicleTrajectory.station)...
            -inputTrajectory(1,fieldsPreSimulationVehicleTrajectory.station)-mdlparam.lookAheadDistance));
        initialAhead.east  = inputTrajectory(initialAhead.index,fieldsPreSimulationVehicleTrajectory.east);
        initialAhead.north = inputTrajectory(initialAhead.index,fieldsPreSimulationVehicleTrajectory.north);
        initial.hdg = atan2(initialAhead.north-initial.north,initialAhead.east-initial.east); %[radians]
        clear initialAhead.index; % clear the variables not in use
        
        initial.speed = inputTrajectory(1,fieldsPreSimulationVehicleTrajectory.speed); %[meter/second]
        initial.wheelspeeds = initial.speed/Wheel_radius*ones(4,1); %[rad/s]
        
        % run the simulation
        fprintf(1,'\nSimulation is running ...\n');

        simlulation_start = tic;
%         profile off
%         profile on -timer 'performance'
        sim('ITSC2021_CoupledPathFollowingSim_discrete.slx', sim_duration_from_aimsun_time);
%       sim('CoupledPathFollowingSim_discrete.slx', sim_duration_from_aimsun_time);
%       sim('CoupledPathFollowingSim_continous.slx', sim_duration_from_aimsun_time);
        
%         profile viewer 
        % print to the command window to see that the simulation works by checking
        % the outputs add a function to check the health
        time_elapsed = toc(simlulation_start);
        
        if flag.verbose
            fprintf(1,'Simulation is done.\n\tWall time to run the simulation: %.5f seconds. \n',time_elapsed);
            fprintf(1,'\tDuration of vehicle trajectory in Aimsun: %.5f seconds. \n',sim_duration_from_aimsun_time);
        end % NOTE: Ends flag.verbose
        
        %% step 9: plot the results
        %just debug few times when inserting the data to database, save time

        if flag.doDebug
            % 1. plot aimsun trajectory, input trajectory of ST Trajectory DB Query, desired trajectory and simulink output trajectory
            fcn_plotTrajectory(inputTrajectory,inputTrajectory,Desired_trajectory,P_CG_ENH,preview_friction);
            
            % 2. plot vehicle on the trajectory generated by simulation
            fcn_plotVehicleOnTrajectory(inputTrajectory,P_CG_ENH(:,1),P_CG_ENH(:,2),P_CG_ENH(:,3),delta,Vehicle_Geometry,tout);
             
            % 3. plot time v.s. velocity and station v.s. velocity
            fcn_plotTimeVelocity(tout,command_velocity,inputTrajectory,Desired_trajectory,Vx,Vy,ahead_safe_velocity);
            
            fcn_plotStationVelocity(command_velocity,inputTrajectory,Desired_trajectory,Vx,Vy,ahead_safe_velocity,P_CG)
            
            % 4. plot lateral look ahead offset error and current lateral error
            Desired_traj = inputTrajectory(:,12:13); % reference trajectory
            Actual_traj = P_CG_ENH(:,1:2); % simulated actual trajectory
            
            lateral_error = fcn_lateralError(Desired_traj,Actual_traj);% actual lateral error
            
            fcn_plotLateralOffsetError(tout,error_lookahead);
            
            
            
            station_cg = [0; cumsum(sqrt(sum(diff(P_CG(:,1:2)).^2,2)))]; % calculate station
            fcn_plotLateralOffsetError_sration(station_cg,error_lookahead)
            
            if flag.speed_preview
                sim_results_with_preview.station = station_cg;
                sim_results_with_preview.error_lookahead = error_lookahead;
                sim_results_with_preview.poistion = P_CG_ENH;
                sim_results_with_preview.command_velocity = command_velocity;
                sim_results_with_preview.Vx = Vx;
                sim_results_with_preview.Vy = Vy;
                sim_results_with_preview.V = sqrt(Vx.^2+Vy.^2);
%                 save('sim_results_with_preview.mat', 'sim_results_with_preview')
            else 
                sim_results_no_preivew.station = station_cg;
                sim_results_no_preivew.error_lookahead = error_lookahead;
                sim_results_no_preivew.poistion = P_CG_ENH;
                sim_results_no_preivew.command_velocity = command_velocity;
                sim_results_no_preivew.Vx = Vx;
                sim_results_no_preivew.Vy = Vy;
                sim_results_no_preivew.V = sqrt(Vx.^2+Vy.^2);
%                 save('sim_results_no_preivew.mat', 'sim_results_no_preivew')
            end
            %% 5. plot lateral look ahead offset error with and without preview data
            % error_lookahead_without_preview = error_lookahead; %lateral error
            % save('error_lookahead_without_preview.mat','error_lookahead_without_preview')
            
            load('sim_results_no_preivew.mat','sim_results_no_preivew')
            load('sim_results_with_preview.mat','sim_results_with_preview')
            fcn_plotLateralOffsetError_both(sim_results_with_preview.error_lookahead,sim_results_no_preivew.error_lookahead,sim_results_no_preivew.station);
            fcn_plotspeed_both(sim_results_with_preview.V,sim_results_no_preivew.V,sim_results_no_preivew.station,sim_results_no_preivew.error_lookahead);
            fcn_plotspeed_cmd_actual(sim_results_with_preview.V,sim_results_with_preview.command_velocity,sim_results_no_preivew.station,sim_results_no_preivew.error_lookahead);
            
            
            %% 6.plot steering angle
            fcn_plotSteeringAngle(delta,tout);
            
            % 7. plot: vehicle CG trajectory and wheels trajectory
            fcn_plotCgAndWheelsTrajectory(frictionGrid_table,lanesCenter_table,P_CG,P_FL,P_FR,P_RL,P_RR)
            
            % 8. plot: vehicle yaw of aimsun trajectory, input trajectory of ST Trajectory DB Query, desired trajectory,  and simulink trajectory 
            fcn_plotStationYaw(queryVehicleTrajectory,inputTrajectory,Desired_trajectory,P_CG)
            
            % 9. plot: vehicle yaw of aimsun trajectory, input trajectory of ST Trajectory DB Query, desired trajectory,  and simulink trajectory 
            fcn_plotTimeYaw(tout,queryVehicleTrajectory,inputTrajectory,Desired_trajectory,P_CG)
            
            % 10. plot: vehicle CG station and friction coefficient at each wheels
            %% fcn_plotStationFrictionCoefficient(P_CG,friction_coeffs,friction_coeffs_estimate) % friction_coeffs is the queried friction
            if flag.noise
                noisy_type = 2; % the estimated friction is noisy 
            else
                noisy_type = 1; % add white Gaussian noise to estimation results
            end
            if noisy_type ==1 % white Gaussian noise
                Signal_noise_ratio = 40; %unit: dB
                friction_coeffs_estimate_noisy = awgn(friction_coeffs_estimate,Signal_noise_ratio,'measured'); % Add white Gaussian noise to signal
            elseif noisy_type ==2 % vehicle sensor 
                friction_coeffs_estimate_noisy = friction_coeffs_estimate;
            end 
            fcn_plotStationFrictionCoefficient(P_CG,friction_coeffs,friction_coeffs_estimate_noisy,'rr') % true and estimated friction without noisy
            fcn_plotStationFrictionCoefficient(P_CG,friction_coeffs,friction_coeffs_estimate_noisy,'fl') % only front left tire
            fcn_plotStationFrictionCoefficient(P_CG,friction_coeffs)% true friction
            
            fcn_plotStationFrictionCoeffEstimateNoisy(P_CG,friction_coeffs_estimate_noisy)
            %%
            error_trueFri_estiFri = mean(friction_coeffs - friction_coeffs_estimate);
            if flag.verbose
                fprintf(1,'The error between true friction and estimated friction for 4 wheels: %.5f,%.5f,%.5f,%.5f. \n',error_trueFri_estiFri);
            end
            % 11. plot: vehicle CG trajectory, wheels trajectory, friction grid and
            % friction coefficient at wheels 3D
            %%fcn_plotVehicleTrajectoryAndFrictionGrid(P_CG,P_FL,P_FR,P_RL,P_RR,lanesCenter_table,frictionGrid_table,friction_coeffs)
            fcn_plotVehicleTrajectoryAndFrictionGrid(P_CG,P_FL,P_FR,P_RL,P_RR,lanesCenter_table,frictionGrid_table,friction_coeffs_estimate)
            
            % 12 plot: station vs preview friction, station vs current/preview friction, time vs current/preview friction
            fcn_plotStationFrictionPreview(P_CG,preview_friction)
            fcn_plotStationCurrentAndPreviewFriction(P_CG,friction_coeffs,preview_friction,'fl')        
            fcn_plotStationCurrentAndPreviewFriction(P_CG,friction_coeffs,preview_friction)   
            
            %fcn_plotTimeCurrentAndPreviewFriction(tout,friction_coeffs,preview_friction,'fl')  
            
            % 13. plot: station vs. preview road curvature
            fcn_plotStationRoadRadius(P_CG,sectionRef_table,road_radius_preview)
            
            % 14. time vs. wheel torques and station vs. wheel torques                
            %wheel_torques_without_preview = wheel_torques;
            % save('wheel_torques_without_preview.mat','wheel_torques_without_preview')
            if flag.controller
                load('wheel_torques_without_preview.mat','wheel_torques_without_preview')
                fcn_plotStationWheelTorque(P_CG,wheel_torques,wheel_torques_without_preview)
            end
            fcn_plotTimeWheelTorque(tout,wheel_torques)
           
        end
       
    end % NOTE: Ends for indexTrajectoryOfInterest
    


return



%% play a music at the end of a code 
load chirp
sound(y,Fs)
load handel
sound(y,Fs)