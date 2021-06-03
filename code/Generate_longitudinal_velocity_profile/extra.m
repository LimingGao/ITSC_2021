clear all
close all
%% try numerical soultion
ft = linspace(0,5,25);
f = ft.^2 - ft - 3;

gt = linspace(1,6,25);
g = 3*sin(gt-0.25);

miu_g = 0.9 * 9.8;

% length = 45;
station_span = [0:0.01:Ls];
% 
Ux_0 = sqrt(miu_g /(slope * Ls) );
% Ux_0 = 20;
f = @(s,Ux)(1/Ux)* sqrt(miu_g^2 - (slope*(Ls-s)*Ux^2)^2);

[s,Ux_frictionCircle] = ode45(f, station_span,Ux_0);
Ux_curceLimit = sqrt(miu_g ./(slope .* (Ls-station_span)) );

Ux_frictionCircle_real = real(Ux_frictionCircle);
ax = (1./Ux_frictionCircle_real).* sqrt((miu_g.^2 - (slope*(Ls-station_span)'.*Ux_frictionCircle_real.^2).^2));
ay = slope*(Ls-station_span)'.*Ux_frictionCircle_real.^2;

% ay = sqrt(miu_g.^2 - ax.^2);
h_fig = figure(210);
set(h_fig,'Name','path');
clf;
hold on
plot(station_span,Ux_curceLimit,'b','LineWidth',2)
plot(station_span,Ux_frictionCircle_real,'g','LineWidth',2)
grid on
xlabel('station')
ylabel('velocity (m/s)')
ylim([20 70])
% axis equal
legend('Curve limit speed','vehicle dynamics limit speed')
box on

h_fig = figure(212);
set(h_fig,'Name','path');
clf;
hold on
% plot(station_span,Ux_curceLimit,'b','LineWidth',2)
plot(station_span,Ux_frictionCircle_real,'g','LineWidth',2)
grid on
xlabel('station')
ylabel('velocity (m/s)')
% ylim([20 70])
% axis equal
legend('Curve limit speed')
box on
%%
h_fig = figure(211);
set(h_fig,'Name','path');
clf;
hold on
plot(station_span,ax,'g','LineWidth',2)
plot(station_span,ay,'b','LineWidth',2)
plot(station_span,sqrt(ay.^2 + ax.^2),'r','LineWidth',2)
grid on
xlabel('station')
ylabel('velocity (m/s)')
legend('ax','ay')
% axis equal
box on


%% 

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
