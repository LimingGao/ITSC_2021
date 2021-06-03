%%%%%%%%%%%%  Script ControllerDesign.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Purpose: Generate the limit speed profile given a path 
% Matlab work Path: ~\GitHub\forgetfulDBs\ControllerDesign

% Author:       Liming
% Created Date: 2021-04-01
%
% Reference: 

%

% To do list:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step 1: Prepare the workspace
% Clear workspace, figures, and console
clear all;  %#ok<CLALL> % Clears workspace to remove any old variables
close all;  % Close all open figures
clc;        % Clear console space

addpath('../Utilities/');    % all the functions and wrapper class
addpath('../DataFiles/');    % all the .mat data files
addpath(genpath('../Figures/'));    % all the image files 

dir.datafiles = ['.' filesep 'DataFiles' filesep]; % all the .mat data files
flag.friction = false; % true: use changing friction 
flag.grade = false; % true: use non-zero grade 
flag.curvature = true; %true: use non-zero curvature 

%%Initialize parameters that commonly change

lamuda = 0.95;
speed_limit = 31.1; %m/s
default_friction = 0.1; % used for constant friction assumption

g = 9.81; %m/s^2

%% step 2: load a path including position, curvature and hypothetical friction and grade
load('path.mat','path');
load('I99_StateCollege73_to_Altoona33_mergedDataNoJumps_20210123.mat','I99_StateCollege73_to_Altoona33_mergedDataNoJumps');

% rename variable 
I99_SC2Altoona = I99_StateCollege73_to_Altoona33_mergedDataNoJumps;
clear I99_StateCollege73_to_Altoona33_mergedDataNoJumps


% show the raw data

h_fig = figure(45);
set(h_fig,'Name','raw path');
clf;
hold on
plot(I99_SC2Altoona.xEast,I99_SC2Altoona.yNorth,'b.','LineWidth',2)
% plot(path.x(1),path.y(1),'r','MarkerSize',20)
grid on
xlabel('xEast')
ylabel('yNorth')
axis equal

%% process the curvature 

% calculate curvature according to the position 
[~, ~, ~, ~,R_spiral,UnitNormalV,concavity]=fnc_parallel_curve(I99_SC2Altoona.xEast, I99_SC2Altoona.yNorth, 1, 0,1,250);

% calcualte the curvature according to the velocity 
%{
yaw_rate = [0; diff(I99_SC2Altoona.Yaw_deg)./diff(I99_SC2Altoona.GPS_Time)];

figure(23)
clf
hold on
% plot(I99_StateCollege73_to_Altoona33_mergedDataNoJumps.station,I99_StateCollege73_to_Altoona33_mergedDataNoJumps.altitude)
plot(I99_SC2Altoona.station,I99_SC2Altoona.Yaw_deg,'b','LineWidth',1)
plot(I99_SC2Altoona.station,yaw_rate,'r','LineWidth',1)

grid on
box on
xlabel('station (m)')
ylabel('yaw and yaw rate (deg)')
% ylim([0 0.01])


figure(24)
clf
hold on
Ux = I99_SC2Altoona.velEast.*cosd(I99_SC2Altoona.Yaw_deg) + I99_SC2Altoona.velNorth.*sind(I99_SC2Altoona.Yaw_deg);
plot(I99_SC2Altoona.GPS_Time,Ux,'b','LineWidth',1)
plot(I99_SC2Altoona.GPS_Time,I99_SC2Altoona.velNorth,'g')
plot(I99_SC2Altoona.GPS_Time,I99_SC2Altoona.velEast,'r','LineWidth',1)
grid on
box on
xlabel('time (s)')
ylabel('velocity (m/s)')
% ylim([0 0.01])

curvature_ss  = (yaw_rate*pi/180)./Ux; % steady state
%}
figure(22)
clf
hold on
% plot(I99_StateCollege73_to_Altoona33_mergedDataNoJumps.station,curvature_ss,'g','LineWidth',1)
plot(I99_SC2Altoona.station,abs(concavity).*1./R_spiral,'r','LineWidth',1)
grid on
box on
xlabel('Station (m)')
ylabel('Curvature')
ylim([-0.01 0.04])



%% choose segment of I99
path_index = find(I99_SC2Altoona.station>10000 & I99_SC2Altoona.station< 50000);

figure(220)
clf
geoplot(I99_SC2Altoona.latitude,I99_SC2Altoona.longitude,'b', ...
    I99_SC2Altoona.latitude(path_index),I99_SC2Altoona.longitude(path_index),'r',...
    'LineWidth',2)

% geolimits([45 62],[-149 -123])
% legend('mergedByKFData.GPS\_Hemisphere','mergedByKFData.MergedGPS','mergedDataNoJumps.MergedGPS')
geobasemap satellite


path.x = I99_SC2Altoona.xEast(path_index);
path.y = I99_SC2Altoona.yNorth(path_index);
path.z = I99_SC2Altoona.zUp(path_index);
path.latitude = I99_SC2Altoona.latitude(path_index);
path.longitude = I99_SC2Altoona.longitude(path_index);
path.curvature = abs(concavity(path_index)).*1./R_spiral(path_index);
path.station = [0; cumsum(sqrt(diff(path.x).^2+diff(path.y).^2))];
% plot path
path_length = max(path.station);

h_fig = figure(45);
set(h_fig,'Name','path');
clf;
hold on
plot(path.x,path.y,'b.','LineWidth',2)
plot(path.x(1),path.y(1),'r','MarkerSize',20)
% scatter3(path.center_x,path.center_y,path.friction,10, path.friction,'.');
grid on
xlabel('xEast')
ylabel('yNorth')
axis equal

h_fig = figure(4545);
set(h_fig,'Name','path curvature');
clf;
hold on
subplot(3,1,1)
plot(path.station,path.curvature,'k','LineWidth',2)
box on
grid on
xlabel('station(m)')
ylabel('curvature')

subplot(3,1,2)
plot(path.station,[0 ; diff(path.curvature)],'b','LineWidth',2)
box on
grid on
xlabel('station(m)')
ylabel('curvature diff')
subplot(3,1,3)
% curvature_vertex = [0; diff([0 ; diff(path.curvature)])];
curvature_vertex = del2(path.curvature);
plot(path.station,curvature_vertex,'b','LineWidth',2)
box on
grid on
xlabel('station(m)')
ylabel('curvature diff diff')

%% assign a friction value 
friction_coefficient  = zeros(length(path.station),1);
% high friction
index_fhigh = ((path.station < 0.2*path_length) | (path.station > 0.85*path_length));
friction_coefficient(index_fhigh)  = 0.9; 
index_f1 = ( path.station <= 0.3750*path_length & path.station >= 0.200*path_length);
friction_coefficient(index_f1)  = 0.3; 
index_f2 = ( path.station <= 0.700*path_length & path.station >= 0.550*path_length);
friction_coefficient(index_f2)  = 0.15; 
index_f3 = ( path.station <= 0.850*path_length & path.station >= 0.700*path_length);
friction_coefficient(index_f3)  = 0.35; 
index_f4 = ( path.station < 0.550 * path_length & path.station > 0.3750*path_length);
friction_coefficient(index_f4)  = 0.42;

if flag.friction
    path.friction = lamuda*friction_coefficient;
else
    friction_coefficient = default_friction *ones(length(friction_coefficient),1);
    path.friction = friction_coefficient;
end

if ~flag.curvature
    path.curvature = zeros(length(friction_coefficient),1);
end

% assign a grade value 
grade  = zeros(length(path.station),1);
if flag.grade
index_grade1 = ((path.station < 0.15*path_length) & (path.station > 0.05*path_length));
grade(index_grade1)  = 6;  % percentage

index_grade2 = ( path.station < 0.50 * path_length & path.station > 0.40*path_length);
grade(index_grade2)  = -6;


end

path.grade = atan(grade./100); % percentage to radii

%% 
%{
%this is not applicable for I99 right now
index_arc = find(path.curvature == 1/200);
index_line = find(path.curvature == 0);
index_spiral = find(path.curvature > 0 & path.curvature < 1/200);

MarkerSize = 5;

h_fig = figure(451);
set(h_fig,'Name','check friction');
x0=300;
y0=400;
width=540;
height=480;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
subplot(3,1,2)
hold on
% plot(path.station,path.friction,'k','LineWidth',1.5)
plot(path.station(index_line),path.friction(index_line),'g.','MarkerSize',MarkerSize)
plot(path.station(index_arc),path.friction(index_arc),'m.','MarkerSize',MarkerSize)
plot(path.station(index_spiral),path.friction(index_spiral),'b.','MarkerSize',MarkerSize)
grid on
xlabel('Station (m)')
ylabel('Friction Coefficient')
xlim([0 max(path.station)])
ylim([0 1])
box on
legend1 = legend('Line','Spiral','Arc');
set(legend1,...
    'Position',[0.303209840313213 0.931259407330584 0.398888891643947 0.035000000635783],...
    'NumColumns',3);

subplot(3,1,1)
hold on
plot(path.station(index_line),path.curvature(index_line),'g.','MarkerSize',MarkerSize)
plot(path.station(index_arc),path.curvature(index_arc),'m.','MarkerSize',MarkerSize)
plot(path.station(index_spiral),path.curvature(index_spiral),'b.','MarkerSize',MarkerSize)
box on
grid on
xlabel('Station(m)')
ylabel('Curvature ')
xlim([0 max(path.station)])
% annotation(h_fig,'textarrow',[0.484444444444445 0.450370370370371],...
%     [0.546666666666667 0.487500000000001],'String',{'spiral'});

subplot(3,1,3)
hold on
% plot(path.station,grade,'k','LineWidth',1.5)
plot(path.station(index_line),grade(index_line),'g.','MarkerSize',MarkerSize)
plot(path.station(index_arc),grade(index_arc),'m.','MarkerSize',MarkerSize)
plot(path.station(index_spiral),grade(index_spiral),'b.','MarkerSize',MarkerSize)
box on
grid on
xlim([0 max(path.station)])
ylim([-10 10])
xlabel('Station(m)')
ylabel('Grade (%) ')



h_fig = figure(45781);
set(h_fig,'Name','check friction');
x0=300;
y0=400;
width=540;
height=480;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
subplot(3,1,1)
hold on
plot(path.station,path.friction,'k','LineWidth',1.5)
grid on
xlabel('Station (m)')
ylabel('Friction Coefficient')
xlim([0 max(path.station)])
ylim([0 1])

subplot(3,1,2)
hold on
plot(path.station,path.curvature,'k','LineWidth',1.5)
box on
grid on
xlabel('Station(m)')
ylabel('Curvature ')
xlim([0 max(path.station)])

subplot(3,1,3)
hold on
plot(path.station,grade,'k','LineWidth',1.5)
box on
grid on
xlim([0 max(path.station)])
ylim([-12 12])
xlabel('Station(m)')
ylabel('Grade (%) ')

h_fig = figure(455);
set(h_fig,'Name','path with friction ');
clf;
hold on
%plot(path.center_x,path.center_y,'b','LineWidth',2)
scatter3(path.x,path.y,path.friction,10, path.friction,'.');
grid on
xlabel('xEast')
ylabel('yNorth')
% axis equal
colormap(jet);
colorbar
box on

%% step2.5 cut the path into segments by curvature, friction, and grade.
if 1==0
ind_curBound =[1; find(abs(curvature_vertex)>1*10^(-6)); length(path.curvature)];
nb_curSeg = length(ind_curBound) -1; % numbers of segments according to the curvature
h_fig = figure(4425);
set(h_fig,'Name','path curvature');
clf;
hold on
% subplot(3,1,1)
plot(path.station,path.curvature,'k.','LineWidth',2)
plot(path.station(ind_curBound),path.curvature(ind_curBound),'r.','MarkerSize',20)
box on
grid on
xlabel('station(m)')
ylabel('curvature')


for i= 1:nb_curSeg
    ind_start = max(ind_curBound(i)-1,1);
    ind_end = ind_curBound(i+1);
    path.segment_num(ind_start:ind_start) = i;
    % p = p1*x + p2
    p = polyfit(path.station(ind_start:ind_end),abs(path.curvature(ind_start:ind_end)),1);  % linear fit the curvature at each segment
    if abs(p(1)) < 10^(-8) && abs(p(2)) < 10^(-8)
        path.segment_type(ind_start:ind_end) = 1; % linear
    elseif abs(p(1)) < 10^(-8) && abs(p(2)) > 10^(-8)
        path.segment_type(ind_start:ind_end) = 2; % circle arc
    elseif p(1) > 10^(-8)
        path.segment_type(ind_start:ind_end) = 3; % spiral slow down
    elseif p(1) < -10^(-8)   
        path.segment_type(ind_start:ind_end) = 4; % spiral speed up
    else 
        error('unknown segment!')
    end
    
end

% for i= 1:nb_segment
%     pathSegment.friction(i,1) = mean(path.friction(boundary_index(i):boundary_index(i+1)-1));
%     pathSegment.U_max(i,1) = min(path.U_max(boundary_index(i):boundary_index(i+1)-1));
%     pathSegment.station_start(i,1) = path.station(boundary_index(i));
% %     boundary_index(i+1)-1
%     pathSegment.station_end(i,1) = path.station(boundary_index(i+1)-1);
%     pathSegment.nb_index(i,1) = boundary_index(i+1) - boundary_index(i);
%     pathSegment.radius(i,1) = 1/mean(path.curvature(boundary_index(i):boundary_index(i+1)-1));
% end

end

%}
%% step 3: plan velocity profile 

% 1. find the curve speed limit of the road given the friction and geomtry

U_max = sqrt(path.friction.*g.*(1./path.curvature));
U_max(U_max>speed_limit) = speed_limit;

path.U_max= U_max;


diff2_curve_limit = [0;0; diff(diff(U_max))];
index_bound = find(abs(diff2_curve_limit)>0.005);

h_fig = figure(452);
set(h_fig,'Name','U max');
clf;
subplot(2,1,1)
hold on
plot(path.station,sqrt(path.friction.*g.*(1./path.curvature)),'b','LineWidth',2)
plot(path.station,U_max,'r','LineWidth',2)
% plot(path.station(index_bound),U_max(index_bound),'k.','MarkerSize',20)

xlim([0 max(path.station)])
ylim([10 speed_limit])
grid on
box on
subplot(2,1,2)
hold on
plot(path.station,[0;0; diff(diff(U_max))],'k','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
box on
xlim([0 max(path.station)])

%% 2. find maximum available speed of vehicle give the vehicle tire dynamic 

delta_s = diff(path.station);

Ux_forwardInt = zeros(1,length(path.station));
Ux_forwardInt(1) = path.U_max(1);
for i = 1:length(delta_s)
    Ux_forwardInt(i+1) = Ux_forwardInt(i) + (delta_s(i)/Ux_forwardInt(i))* (sqrt((path.friction(i)*g*cos(path.grade(i)))^2 - (path.curvature(i)*Ux_forwardInt(i)^2)^2) - g*sin(path.grade(i)));
    if Ux_forwardInt(i+1)>U_max(i+1)
       Ux_forwardInt(i+1)=U_max(i+1);
    end
end


Ux_backwardInt = zeros(1,length(path.station));
Ux_backwardInt(end) = path.U_max(end);

for i = length(delta_s):-1:1
    Ux_backwardInt(i) = Ux_backwardInt(i+1) + (delta_s(i)/Ux_backwardInt(i+1))* (sqrt((path.friction(i+1)*g*cos(path.grade(i+1)))^2 - (path.curvature(i+1)*Ux_backwardInt(i+1)^2)^2) + g*sin(path.grade(i+1)));
%     if Ux_backwardInt(i)>Ux_forwardInt(i)
%         Ux_backwardInt(i)=Ux_forwardInt(i);
%     end
    
    if Ux_backwardInt(i)>U_max(i)
        Ux_backwardInt(i)=U_max(i);
    end
end

v = min(real(Ux_forwardInt), real(Ux_backwardInt));
ax = [0 0.5.*diff(v.^2)./diff(path.station')];
ay = v.^2.*path.curvature';
a_max = cos(path.grade').*g.*path.friction'./lamuda + g.*abs(sin(path.grade'));

h_fig = figure(4520);
set(h_fig,'Name','U max');
x0=300;
y0=400;
width=540;
height=300;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on

% plot(path.station,sqrt(path.friction.*g.*(1./path.curvature)),'b','LineWidth',2)
plot(path.station,U_max,'b','LineWidth',2)
plot(path.station,Ux_forwardInt,'--','Color',[0 1 1],'LineWidth',2.5)
plot(path.station,Ux_backwardInt,'m-.','LineWidth',2.5)
plot(path.station,min(Ux_forwardInt, Ux_backwardInt),'g-','LineWidth',1.5)
% plot(755.6,23.63,'r.','MarkerSize',25) % C
plot(900.4,50,'k.','MarkerSize',25) % C'
plot(900.4,38.,'r.','MarkerSize',25) % C
plot(1345.5,16.713,'r.','MarkerSize',25 ) % D
plot(419,23.63,'r.','MarkerSize',25)% B
plot(315,40.9,'r.','MarkerSize',25)% A
legend1 = legend('curve limit speed','forward limit','backward limit', 'combined results ');
set(legend1,...
    'Position',[0.130432058173258 0.145111107296414 0.540000008371141 0.121666669209798],...
    'NumColumns',2);
grid on
box on
xlabel('Station (m)')
ylabel('Ux (m/s)')
xlim([0 max(path.station)])
ylim([5 speed_limit+5])
% xlim([0 499])
%% safe and stable speed 
path.velocity_unsafe = 25*ones(size(path.station));
path.velocity_safe = min(path.velocity_unsafe,v') ;
path.friction_true = friction_coefficient;
filename = 'path_desired.mat';
% save(filename,'path')
% lower_bound = v;
% save('lower_bound.mat','lower_bound')
load('lower_bound.mat')


h_fig = figure(455077);
set(h_fig,'Name','U_profile_second');
x0=300;
y0=400;
width=540;
height=300;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b--','LineWidth',2)
plot(path.station,v,'g','LineWidth',2)
% plot(path.station(index_plan),v(index_plan),'r.','LineWidth',2)
% plot(path.station(index_plan_jump),v(index_plan_jump),'k.','MarkerSize',20)
plot(path.station,path.velocity_safe,'m--','LineWidth',2)
% plot(path.station,lower_bound,'r-','LineWidth',2)


grid on
xlabel('Station (m)')
ylabel('Ux (m/s)')
xlim([0 max(path.station)])
ylim([0 speed_limit+5])
legend('Curve limit speed','Speed profile for actual ¦Ì','Speed profile for safe')

% legend('Curve limit speed','Speed profile for actual ¦Ì','Speed profile for confident low bound ¦Ì = 0.2')
box on
%% calcualet the preview distance 
% v = min(real(Ux_forwardInt), real(Ux_backwardInt));
% difference_max_profile = U_max - v';  % thid may not be the correct preview distance 
Ux_backwardInt = real(Ux_backwardInt);
difference_max_profile = U_max - Ux_backwardInt';
index_plan = find(difference_max_profile>0.1);

if ~isempty(index_plan)
    
    diff_index_plan = diff(index_plan);
    diff_index_plan_jump = [1 ;find(abs(diff_index_plan)>1) ;find(abs(diff_index_plan)>1)+1; length(diff_index_plan)];
    index_plan_jump = sort(index_plan(diff_index_plan_jump));
    
    
    d_preview= zeros(size(path.station));
    
    for i = 1:1:(length(index_plan_jump)/2)
        for j = 1:1:floor((index_plan_jump(2*i)-index_plan_jump(2*i-1)))
            d_preview(index_plan_jump(2*i-1)+j) = path.station((index_plan_jump(2*i))) - path.station((index_plan_jump(2*i-1)+j));
            %d_preview(index_plan_jump(i)+j) = path.station((index_plan_jump(i+1))- path.station((index_plan_jump(i)+j)));
        end
    end
else
    index_plan_jump = [];
    d_preview= zeros(size(path.station));
end
h_fig = figure(4550);
set(h_fig,'Name','U_profile_second');
x0=300;
y0=400;
width=540;
height=250;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b--','LineWidth',2)
plot(path.station,v,'g','LineWidth',2)
plot(path.station(index_plan),Ux_backwardInt(index_plan),'r.','LineWidth',2)
plot(path.station(index_plan_jump),Ux_backwardInt(index_plan_jump),'k.','MarkerSize',20)
plot(path.station,Ux_backwardInt,'m--','LineWidth',2)
plot(200,35,'m.','MarkerSize',20)

grid on
xlabel('Station (m)')
ylabel('Ux (m/s)')
xlim([0 max(path.station)])
ylim([0 speed_limit+5])
legend('Curve limit speed','Speed profile for actual ¦Ì','Speed profile for safe')

box on
%
h_fig = figure(8078);
set(h_fig,'Name','diff');
x0=300;
y0=400;
width=540;
height=250;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
plot(path.station,d_preview,'k','LineWidth',2)
grid on
xlabel('Station (m)')
ylabel('Preview Distance (m)')
xlim([0 max(path.station)])

box on
%%
index_dangerous = find(d_preview>50);

figure(220)
clf
% geoplot(I99_SC2Altoona.latitude,I99_SC2Altoona.longitude,'b', ...
%     I99_SC2Altoona.latitude(path_index),I99_SC2Altoona.longitude(path_index),'r',...
%     'LineWidth',2)
geoplot(I99_SC2Altoona.latitude(path_index),I99_SC2Altoona.longitude(path_index),'b', ...
    path.latitude(index_dangerous),path.longitude(index_dangerous),'r.',...
    'LineWidth',2,'MarkerSize',25)
% geolimits([45 62],[-149 -123])
legend('path','path dangerous region')
geobasemap satellite

figure(2222)
clf
hold on
% plot(I99_StateCollege73_to_Altoona33_mergedDataNoJumps.station,curvature_ss,'g','LineWidth',1)
plot(path.station,path.curvature,'b','LineWidth',1)
plot(path.station(index_dangerous),path.curvature(index_dangerous),'r.','LineWidth',1)
grid on
legend('all curvature','curvature for ¦Ì = 0.1 and d_p > 50m')

box on
xlabel('Station (m)')
ylabel('Curvature')
ylim([-0.001 0.004])
