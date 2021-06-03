%%%%%%%%%%%%  Script ControllerDesign.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Purpose: Vehicle Dynamics Analysis and Controller Design
% Matlab work Path: ~\GitHub\forgetfulDBs\ControllerDesign

% Author:       Liming
% Created Date: 2020-05-15
%
% Reference: 
% [1] Craig E. Beal & Christina Boyd (2019) Coupled
%     lateral-longitudinal vehicle dynamics and control design with
%     three-dimensional state portraits, Revisions: 2020-05-27:
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


%%Initialize parameters that commonly change
U = 30; % m/s
steering_amplitude = 20;
%%Initialize constant vehicle parameters
vehicle(1).m          = 1031.92; % kg
vehicle(1).Iz         = 1850.5; % kg-m^2
vehicle(1).a          = 0.9271; % meters
vehicle(1).b          = 1.5621;  % meters
vehicle(1).Caf        = -77500; % N/rad;
%vehicle(1).Car	      = -43000; % N/rad;  
vehicle(1).Car	      = -116250; % N/rad;
g = 9.8; %m/s^2

% used at Craig E. Beal & Christina Boyd (2019) Coupled
% lateral-longitudinal vehicle dynamics and control design with
% three-dimensional state portraits, Vehicle System Dynamics, 57:2,
% 286-313.

vehicle(2).m          = 1640; % kg
vehicle(2).Iz         = 3500; % kg-m^2
vehicle(2).a          = 1.39; % meters
vehicle(2).b          = 1.45;  % meters
vehicle(2).d          = 0.81;  % meters, Half track width
vehicle(2).Cxf        = 40000; % N/rad;
vehicle(2).Caf        = 55000; % N/rad;
vehicle(2).Cxr	      = 50000; % N/rad;
vehicle(2).Car	      = 55000; % N/rad;


%% step 2: load a path preview including position, curvature and 
load('lanesCenter_table_section201.mat','lanesCenter_table');

path = lanesCenter_table(lanesCenter_table.road_lane_id == 20101,:);

h_fig = figure(45);
set(h_fig,'Name','path');
clf;
hold on
plot(path.center_x,path.center_y,'b','LineWidth',2)
% scatter3(path.center_x,path.center_y,path.friction,10, path.friction,'.');
grid on
xlabel('xEast')
ylabel('yNorth')
axis equal

% assign a friction value 
friction_coefficient  = zeros(height(path),1);
% high friction
index_fhigh = ((path.station < 200) | (path.station > 850)|( path.station < 550 & path.station > 350));
friction_coefficient(index_fhigh)  = 0.9; 
index_f1 = ( path.station <= 350 & path.station >= 200);
friction_coefficient(index_f1)  = 0.2; 
index_f2 = ( path.station <= 700 & path.station >= 550);
friction_coefficient(index_f2)  = 0.6; 
index_f3 = ( path.station <= 850 & path.station >= 700);
friction_coefficient(index_f3)  = 0.5; 

path.friction = friction_coefficient;

h_fig = figure(451);
set(h_fig,'Name','check friction');
clf;
hold on
plot(path.station,path.friction,'b.','LineWidth',2)

grid on
xlabel('station')
ylabel('friction')
ylim([0 1])

h_fig = figure(455);
set(h_fig,'Name','path');
clf;
hold on
%plot(path.center_x,path.center_y,'b','LineWidth',2)
scatter3(path.center_x,path.center_y,path.friction,10, path.friction,'.');
grid on
xlabel('xEast')
ylabel('yNorth')
axis equal
colormap(jet);
colorbar

%% step 3: plan to velocity profile 

% 1. find the maximum allowable speed of the road given the friction and geomtry
% 
speed_limit = 35; 
U_max = sqrt(path.friction.*g.*(1./path.curvature));
U_max(U_max>speed_limit) = speed_limit;
% 2. find maximum available speed of vehicle give the vehicle tire dynamic 




path.U_profile= U_max;
h_fig = figure(452);
set(h_fig,'Name','U max');
clf;
hold on
%plot(path.station,path.U_profile,'b','LineWidth',2)
% plot(path.station,1./path.curvature,'m.','LineWidth',2)
plot(station_c,sqrt(friction.*g.*200),'b','LineWidth',2)
[mi,index] = min(abs(Ux-28));
station_c(index)
plot(station_c(1:index)+200-station_c(index),Ux(1:index),'g','LineWidth',2)
plot(station_c+300,Ux2,'g','LineWidth',2)
plot([200,300],[28 28],'g','LineWidth',2)
plot([0,200-station_c(index)],[42 42],'g','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
ylim([20 45])
xlim([0 499])


%% plot the solution of longitudinal velocity
syms  y(t) t a b
eqn = diff(y,t) == (1/y)*sqrt(a^2-y^4/b^2);
S = dsolve(eqn)



station_c = 0:0.1:500;
R = 200;
miu_g = 9.81*0.4;
u_x_initial = 25;

c1 = (R/2)* atan(sqrt(u_x_initial^4/(miu_g^2*R^2 -u_x_initial^4 )))

u_x_acce = sqrt(miu_g*R) .* sqrt(tan(2*(c1+station_c)./R))./((tan(2.*(c1+station_c)./R).^2+1).^(1/4));
% the period of this equation is pi*R/2

u_x_dece = sqrt(miu_g*R) .* sqrt(-tan(2*(-c1+station_c)./R))./((tan(2.*(-c1+station_c)./R).^2+1).^(1/4));


h_fig = figure(452);
set(h_fig,'Name','U max');
clf;
hold on

u_x_acce(real(u_x_acce)<1) = max(real(u_x_acce));
u_x_dece(real(u_x_dece)>25) = min(real(u_x_dece));

plot(station_c,u_x_acce,'b','LineWidth',2)

plot(station_c,u_x_dece,'g','LineWidth',2)
% plot(station,Ux2,'r--','LineWidth',2)
grid on
box on
xlabel('station')
ylabel('u x')
xlim([0 200])
ylim([0 45])
legend('ux\_acce','ux\_dece','ux\_dece\_simulation')

% station from velocity 

c_s = R * atan(R*sqrt(miu_g^2 -u_x_initial^4/R^2)/u_x_initial^2)/2

u_x_final = 15;
if u_x_final> u_x_final
    S_acce = c_s - R * atan(R*sqrt(miu_g^2 -u_x_final^4/R^2)/u_x_final^2)/2
else
    S_dece = -c_s + R * atan(R*sqrt(miu_g^2 -u_x_final^4/R^2)/u_x_final^2)/2
end
%%


c1 = R*sqrt(miu_g^2*R^2-u_x_initial^2);
c2 = -c1;
u_x_dece = sqrt(miu_g^2*R^4 - (station_c + c1).^2)./R;

u_x_acc = sqrt(miu_g^2*R^4 - (-station_c + c1).^2)./R;



%% play a music at the end of a code
% load chirp
% sound(y,Fs)
% load handel
% sound(y,Fs)

% lineHandles = find_system(gcs,'FindAll','On','SearchDepth',1,'Type','Line');
% Simulink.BlockDiagram.routeLine(lineHandles);