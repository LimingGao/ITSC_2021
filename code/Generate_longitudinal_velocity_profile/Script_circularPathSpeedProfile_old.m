%%%%%%%%%%%%  Script Script_circularPathSpeedProfile.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Purpose: 
%       Calculate the speed profile of circular arc  

% Author:       Liming
% Created Date: 2020-05-15
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


g = 9.81; %m/s^2

%% step 2: load or design a path preview including position, curvature and friction 
% load('lanesCenter_table_section201.mat','lanesCenter_table');

R = 200; % Radius 
t= 0:0.001:1.5*pi; % arc angle
% % path.center_x = R*sin(t);
path.center_y = R*cos(t);
path.station = [0 cumsum(sqrt(diff(path.center_x).^2+diff(path.center_y).^2))];
path_station = max(path.station);

[~, ~, ~, ~,R,UnitNormalV]=fnc_parallel_curve(path.center_x, path.center_y, 1, 0,0,0);
path.curvature = 1./R;
h_fig = figure(45);
set(h_fig,'Name','path');
clf;
hold on
plot(path.center_x,path.center_y,'b','LineWidth',2)
plot(path.center_x(1),path.center_y(1),'r','MarkerSize',20)
% scatter3(path.center_x,path.center_y,path.friction,10, path.friction,'.');
grid on
xlabel('xEast')
ylabel('yNorth')
axis equal

% assign a friction value 
friction_coefficient  = zeros(length(path.station),1);
% high friction
index_fhigh = ((path.station < 0.2*path_station) | (path.station > 0.85*path_station));
friction_coefficient(index_fhigh)  = 0.9; 
index_f1 = ( path.station <= 0.3750*path_station & path.station >= 0.200*path_station);
friction_coefficient(index_f1)  = 0.3; 
index_f2 = ( path.station <= 0.700*path_station & path.station >= 0.550*path_station);
friction_coefficient(index_f2)  = 0.2; 
index_f3 = ( path.station <= 0.850*path_station & path.station >= 0.700*path_station);
friction_coefficient(index_f3)  = 0.4; 
index_f4 = ( path.station < 0.550 * path_station & path.station > 0.3750*path_station);
friction_coefficient(index_f4)  = 0.55; 
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
box on
%% step2.5 cut the path into segments by friction vlaues 
[counts,edges] = histcounts(path.friction,100,'BinLimits',[0,1]); % Calculate a 101-bin histogram for the friction.
h_fig = figure(15565);
set(h_fig,'Name','histogram of friction data');
histogram(path.friction,edges)
%stem(edges,counts)
% execute the k-means cluster
nbs_cluster = length(counts(counts>10));
[cluster_idx,cluster_centroid_locations,sumDist] = kmeans(path.friction,nbs_cluster,'Replicates',1,'Display','final');

cluster_idx_diff = [0 ;diff(cluster_idx)];
[row,~,v] = find(cluster_idx_diff~=0);
boundary_index = [1 ;row;length(path.station)+1];
nb_segment = length(boundary_index) -1;
path.segment = zeros(length(path.station),1);
for i= 1:nb_segment
    path.segment(boundary_index(i):boundary_index(i+1)) = i;
end

%% step 3: plan velocity profile 

% 1. find the maximum allowable speed of the road given the friction and geomtry
% 
speed_limit = 35; 
U_max = sqrt(path.friction.*g.*(1./path.curvature));
U_max(U_max>speed_limit) = speed_limit;
% 2. find maximum available speed of vehicle give the vehicle tire dynamic 

path.U_max= U_max;
h_fig = figure(452);
set(h_fig,'Name','U max');
clf;
hold on
%plot(path.station,path.U_profile,'b','LineWidth',2)
% plot(path.station,1./path.curvature,'m.','LineWidth',2)
plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,U_max,'r','LineWidth',2)
% [mi,index] = min(abs(Ux-28));
% station_c(index)
% plot(station_c(1:index)+200-station_c(index),Ux(1:index),'g','LineWidth',2)
% plot(station_c+300,Ux2,'g','LineWidth',2)
% plot([200,300],[28 28],'g','LineWidth',2)
% plot([0,200-station_c(index)],[42 42],'g','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
ylim([10 45])
% xlim([0 499])


%% step4: using friction ellipse 
%find minimum velocity for each segment
for i= 1:nb_segment
    pathSegment.friction(i,1) = mean(path.friction(boundary_index(i):boundary_index(i+1)-1));
    pathSegment.U_max(i,1) = min(path.U_max(boundary_index(i):boundary_index(i+1)-1));
    pathSegment.station_start(i,1) = path.station(boundary_index(i));
%     boundary_index(i+1)-1
    pathSegment.station_end(i,1) = path.station(boundary_index(i+1)-1);
    pathSegment.nb_index(i,1) = boundary_index(i+1) - boundary_index(i);
    pathSegment.radius(i,1) = 1/mean(path.curvature(boundary_index(i):boundary_index(i+1)-1));
end

% find the segments which need to slow down
reduceSpeed_idx = find([diff(pathSegment.U_max)<0]>0);

% plan the velocity profile for segments which need to slow down ()
path.U_profile_intial= path.U_max;
for i = 1:length(reduceSpeed_idx)
    u_x_initial = pathSegment.U_max(reduceSpeed_idx(i)+1);
    miu_g = pathSegment.friction(reduceSpeed_idx(i))*g;
    R = pathSegment.radius(reduceSpeed_idx(i));
    
    c1 = (R/2)* atan(sqrt(u_x_initial^4/(miu_g^2*R^2 -u_x_initial^4 )));
    
    station_c = linspace(pathSegment.station_start(reduceSpeed_idx(i)),pathSegment.station_end(reduceSpeed_idx(i)),pathSegment.nb_index(reduceSpeed_idx(i))) - pathSegment.station_start(reduceSpeed_idx(i));
    u_x_acce = sqrt(miu_g*R) .* sqrt(tan(2*(c1+station_c)./R))./((tan(2.*(c1+station_c)./R).^2+1).^(1/4));
    u_x_acce(real(u_x_acce)<1) = max(real(u_x_acce));
    
    % Umax limitation
    u_x_acce(u_x_acce>pathSegment.U_max(reduceSpeed_idx(i))) = pathSegment.U_max(reduceSpeed_idx(i));
    
    % plot
    h_fig = figure(2);
    set(h_fig,'Name','U max');
    clf;
    hold on
    plot(station_c,u_x_acce,'b','LineWidth',2)
    grid on
    
    path.U_profile_intial(boundary_index(reduceSpeed_idx(i)):boundary_index(reduceSpeed_idx(i)+1)-1) = flip(u_x_acce);
end


h_fig = figure(450);
set(h_fig,'Name','U max');
clf;
hold on

% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
ylim([10 45])
legend('Curve limit speed')
box on


h_fig = figure(453);
set(h_fig,'Name','U max');
clf;
hold on

% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b','LineWidth',2)
plot(path.station,path.U_profile_intial,'c','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
ylim([10 45])
legend('Curve limit speed','Speed profile considering deceleration')
box on


%% (TO DO LIST)check no segment need reduce speed

% find the segments which need to speed up
SpeedUp_idx = find([0;diff(pathSegment.U_max)>0]>0);
path.U_profile_second= path.U_profile_intial;
for i = 1:length(SpeedUp_idx)
    u_x_initial = pathSegment.U_max(SpeedUp_idx(i)-1);
    miu_g = pathSegment.friction(SpeedUp_idx(i))*g;
    R = pathSegment.radius(SpeedUp_idx(i));
    
    c1 = (R/2)* atan(sqrt(u_x_initial^4/(miu_g^2*R^2 -u_x_initial^4 )));
    
    station_c = linspace(pathSegment.station_start(SpeedUp_idx(i)),pathSegment.station_end(SpeedUp_idx(i)),pathSegment.nb_index(SpeedUp_idx(i))) - pathSegment.station_start(SpeedUp_idx(i));
    u_x_acce = sqrt(miu_g*R) .* sqrt(tan(2*(c1+station_c)./R))./((tan(2.*(c1+station_c)./R).^2+1).^(1/4));
    u_x_acce(real(u_x_acce)<1) = max(real(u_x_acce));
    
    % Umax limitation
    u_x_acce(u_x_acce>pathSegment.U_max(SpeedUp_idx(i))) = pathSegment.U_max(SpeedUp_idx(i));
    
    % plot
    h_fig = figure(3);
    set(h_fig,'Name','U max');
    clf;
    hold on
    plot(station_c,u_x_acce,'b','LineWidth',2)
    grid on
    
    path.U_profile_second(boundary_index(SpeedUp_idx(i)):boundary_index(SpeedUp_idx(i)+1)-1) = u_x_acce;
end

% try backword and forward intergration

delta_s = diff(path.station);

Ux_forwardInt = zeros(1,length(path.station));
Ux_forwardInt(1) = path.U_max(1);
for i = 1:length(delta_s)
    Ux_forwardInt(i+1) = Ux_forwardInt(i) + (delta_s(i)/Ux_forwardInt(i))* sqrt((path.friction(i)*g)^2 - (path.curvature(i)*Ux_forwardInt(i)^2)^2);
    if Ux_forwardInt(i+1)>U_max(i+1)
       Ux_forwardInt(i+1)=U_max(i+1);
    end
end

Ux_backwardInt = zeros(1,length(path.station));
Ux_backwardInt(end) = path.U_max(end);

for i = length(delta_s):-1:1
    Ux_backwardInt(i) = Ux_backwardInt(i+1) + (delta_s(i)/Ux_backwardInt(i+1))* sqrt((path.friction(i+1)*g)^2 - (path.curvature(i+1)*Ux_backwardInt(i+1)^2)^2);
    if Ux_backwardInt(i)>U_max(i)
       Ux_backwardInt(i)=U_max(i);
    end
end



h_fig = figure(4540);
set(h_fig,'Name','U_profile_second');
clf;
hold on
% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b','LineWidth',2)
% plot(path.station,path.U_profile_intial,'c--','LineWidth',2)
plot(path.station,path.U_profile_second,'m-','LineWidth',2)
plot(path.station,Ux_forwardInt,'c--','LineWidth',2)

grid on
box on
xlabel('station')
ylabel('u x')
ylim([10 45])
legend('Curve limit speed','Speed profile considering acceleration')


h_fig = figure(454);
set(h_fig,'Name','U_profile_second');
clf;
hold on

% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b','LineWidth',2)
plot(path.station,Ux_forwardInt,'k-','LineWidth',1.5)
plot(path.station,Ux_backwardInt,'k-','LineWidth',1.5)
plot(path.station,path.U_profile_intial,'c--','LineWidth',2)
plot(path.station,path.U_profile_second,'m--','LineWidth',2)

grid on
xlabel('station')
ylabel('u x')
ylim([10 45])
legend('Curve limit speed','Speed profile considering deceleration','Speed profile considering acceleration')

%% Final Speed profile considering deceleration and acceleration

path.U_profile = min(path.U_profile_intial, path.U_profile_second);

h_fig = figure(4550);
set(h_fig,'Name','U_profile_second');
clf;
hold on

% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
plot(path.station,path.U_max,'b','LineWidth',2)
% plot(path.station,path.U_profile_intial,'c--','LineWidth',2)
% plot(path.station,path.U_profile_second,'m--','LineWidth',2)
plot(path.station,path.U_profile,'g','LineWidth',2)
plot(path.station,min(Ux_forwardInt, Ux_backwardInt),'k--','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
ylim([10 45])
legend('Curve limit speed','Speed profile considering deceleration and acceleration')
box on

h_fig = figure(4590);
set(h_fig,'Name','U_profile_second');
clf;
hold on

% plot(path.station,sqrt(path.friction.*g.*200),'b','LineWidth',2)
% plot(path.station,path.U_max,'b','LineWidth',2)
% plot(path.station,path.U_profile_intial,'c--','LineWidth',2)
% plot(path.station,path.U_profile_second,'m--','LineWidth',2)
plot(path.station,path.U_profile,'g','LineWidth',2)
plot(path.station,15 *ones(size(path.station)),'k--','LineWidth',2)
grid on
xlabel('station')
ylabel('u x')
ylim([10 45])
legend('Speed profile at actual friction','Speed profile at lower friction bound')
box on


%% plot the analytical solution of longitudinal velocity

% syms  y(t) t a b
% eqn = diff(y,t) == (1/y)*sqrt(a^2-y^4/b^2);
% S = dsolve(eqn)

station_c = 0:0.1:500;
R = 200;
miu_g = 9.81*0.8;
u_x_initial = 15;

%%approximation solution 

syms y(x) a b x c
assume(a >= 0)
assume(b >= 0)
f = (1/y)*sqrt(a^2-(b*y^2)^2);
% G = functionalDerivative(f,y)
y_pp = diff(f,x);

y_ppp = diff(y_pp,x);

y_pppp = diff(y_ppp,x);
y_5p = diff(y_pppp,x);

y_6p = diff(y_5p,x);
% argnames(y_pp)

% x = 0
y_p_zero = subs(f,[x y],[0 c]);

y_pp_imd = subs(y_pp,[x diff(y(x), x)],[0 y_p_zero]);
y_pp_zero = subs(y_pp_imd,y(0),c);
simplify(y_pp_zero);

y_ppp_imd = subs(y_ppp,[x diff(y(x), x) diff(y(x), x, x)],[0 y_p_zero y_pp_zero]);
y_ppp_zero = subs(y_ppp_imd,y(0),c);
simplify(y_ppp_zero);


y_pppp_imd = subs(y_pppp,[x diff(y(x), x) diff(y(x), x, x) diff(y(x), x, x,x)],[0 y_p_zero y_pp_zero y_ppp_zero]);
y_pppp_zero = subs(y_pppp_imd,y(0),c);
simplify(y_pppp_zero);

y_5p_imd = subs(y_5p,[x diff(y(x), x) diff(y(x), x, x) diff(y(x), x, x,x) diff(y(x), x, x,x,x)],[0 y_p_zero y_pp_zero y_ppp_zero y_pppp_zero]);
y_5p_zero = subs(y_5p_imd,y(0),c);
simplify(y_5p_zero);

y_6p_imd = subs(y_6p,[x diff(y(x), x) diff(y(x), x, x) diff(y(x), x, x,x) diff(y(x), x, x,x,x) diff(y(x), x, x,x,x,x)],[0 y_p_zero y_pp_zero y_ppp_zero y_pppp_zero y_5p_zero   ]);
y_6p_zero = subs(y_6p_imd,y(0),c);
simplify(y_6p_zero);

% coefficient = [c y_p_zero y_pp_zero/2  y_ppp_zero/6  y_pppp_zero/24];
coefficient = [c y_p_zero y_pp_zero/2  y_ppp_zero/6  y_pppp_zero/24 y_5p_zero/factorial(5) y_6p_zero/factorial(6)];


coefficient_num = double(subs(coefficient,[a b c],[miu_g 1/R u_x_initial]));

% Ux_approximation = coefficient_num(1) + coefficient_num(2).*station_c + coefficient_num(3).*station_c.^2 + ...
%                    coefficient_num(4) .*station_c.^3 + coefficient_num(5) .*station_c.^4;

Ux_approximation = coefficient_num(1) + coefficient_num(2).*station_c + coefficient_num(3).*station_c.^2 + ...
                   coefficient_num(4) .*station_c.^3 + coefficient_num(5) .*station_c.^4 + ...
                   0*coefficient_num(6) .*station_c.^5 +0* coefficient_num(7) .*station_c.^6;


%% closed form solution 


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
plot(station_c,Ux_approximation,'m--','LineWidth',2)

% plot(station,Ux2,'r--','LineWidth',2)
grid on
box on
xlabel('station')
ylabel('u x')
xlim([0 200])
ylim([0 45])
legend('ux\_acce','ux\_dece','ux\_dece\_simulation')

%% distance to acce or dece 


c_s = R * atan(R*sqrt(miu_g^2 -u_x_initial^4/R^2)/u_x_initial^2)/2

u_x_final = 28;
if u_x_final> u_x_initial
    S_acce = c_s - R * atan(R*sqrt(miu_g^2 -u_x_final^4/R^2)/u_x_final^2)/2
else
    S_dece = -c_s + R * atan(R*sqrt(miu_g^2 -u_x_final^4/R^2)/u_x_final^2)/2
end

c_s2 = (R/2)* atan(sqrt(u_x_initial^4/(miu_g^2*R^2 -u_x_initial^4 )))

if u_x_final> u_x_initial
    S_acce_2 = -c_s2 + (R/2)* atan(sqrt(u_x_final^4/(miu_g^2*R^2 -u_x_final^4 )))
else
    S_dece_2 = c_s2 - (R/2)* atan(sqrt(u_x_final^4/(miu_g^2*R^2 -u_x_final^4 )))
end



%% the dece distance with maximum initial velocity 
u_x_final = 0.0000000001;
miu_g_X = 9.81*(0.1:0.05:1);
R_Y = 200:1:201;
[miu_g,R] = meshgrid(miu_g_X,R_Y);

S_dece = R.* atan(R.*sqrt(miu_g.^2 -u_x_final.^4./R.^2)./u_x_final^2)./2;
h_fig = figure(452111);
set(h_fig,'Name','stop distance ');
clf;
hold on
mesh(miu_g,R,S_dece)
colorbar
% plot(miu_g,S_dece,'b','LineWidth',2)

% plot(station_c,u_x_dece,'g','LineWidth',2)
grid on
box on
xlabel('miu\_g')
ylabel('radius')
% xlim([0 200])
% ylim([0 45])
% legend('ux\_acce','ux\_dece','ux\_dece\_simulation')

%% relation between inital speed and stop distance 
R = 200;
miu_g = 0.7* 9.81;
u_x_initial = 1:1: sqrt(miu_g*R);
u_x_final = 1;
c_s = R.* atan(R.*sqrt(miu_g.^2 -u_x_initial.^4./R.^2)./u_x_initial.^2)./2;

S_dece = -c_s + R.* atan(R.*sqrt(miu_g.^2 -u_x_final.^4./R.^2)./u_x_final^2)./2;
h_fig = figure(452191);
set(h_fig,'Name','stop distance ');
clf;
hold on
plot(u_x_initial,S_dece,'b','LineWidth',2)
grid on
box on
xlabel('u\_x\_initial [m/s]')
ylabel('stop distance [m]')

%% relation between radius and stop distance 
R = 200:1:10000;
miu_g = 0.7 * 9.81;
u_x_initial = 37;
u_x_final = 1;
c_s = R.* atan(R.*sqrt(miu_g.^2 -u_x_initial.^4./R.^2)./u_x_initial.^2)./2;

S_dece = -c_s + R.* atan(R.*sqrt(miu_g.^2 -u_x_final.^4./R.^2)./u_x_final^2)./2;
h_fig = figure(4191);
set(h_fig,'Name','stop distance ');
clf;
hold on
plot(R ,S_dece,'b','LineWidth',2)
grid on
box on
xlabel('R  [m]')
ylabel('stop distance [m]')

%% relation between friction and stop distance 
R = 200;
miu_g = (0.2:0.1:1) * 9.81;
u_x_initial = 19.8;
u_x_final = 0.00001;
c_s = R.* atan(R.*sqrt(miu_g.^2 -u_x_initial.^4./R.^2)./u_x_initial.^2)./2;

S_dece = -c_s + R.* atan(R.*sqrt(miu_g.^2 -u_x_final.^4./R.^2)./u_x_final^2)./2;
h_fig = figure(4192);
set(h_fig,'Name','stop distance ');
clf;
hold on
plot(miu_g/9.81 ,S_dece,'b','LineWidth',2)
grid on
box on
xlabel('friction coeficient [m]')
ylabel('stop distance [m]')


%%


c1 = R*sqrt(miu_g^2*R^2-u_x_initial^2);
c2 = -c1;
u_x_dece = sqrt(miu_g^2*R^4 - (station_c + c1).^2)./R;

u_x_acc = sqrt(miu_g^2*R^4 - (-station_c + c1).^2)./R;


%%
Rc = 300;

miu_conf_g = (0.2)*9.81;

for i=1:length(miu_conf_g)
u_x_ini = 0:1:45;

La = Rc*pi/4-(Rc/2)* atan(miu_conf_g(i)*Rc./sqrt(u_x_ini.^4 - (miu_conf_g(i)*Rc).^2 ));

index = find(u_x_ini<sqrt(miu_conf_g(i)*Rc));

La(index) = 0;

h_fig = figure(4192);
set(h_fig,'Name','stop distance ');
% clf;
hold on
plot(u_x_ini,La,'-','LineWidth',1.5)
grid on
box on
xlabel('current speed [m/s]')
ylabel('Preview Distance[m]')

end

legend('¦Ì_c_o_n_f = 0.2','¦Ì_c_o_n_f = 0.4','¦Ì_c_o_n_f = 0.6','¦Ì_c_o_n_f = 0.8')


h_fig = figure(418892);
set(h_fig,'Name','stop distance ');
% clf;
hold on
plot(u_x_ini,La,'-','LineWidth',1.5)
grid on
box on
xlabel('current speed [m/s]')
ylabel('Preview Distance[m]')

legend('R_c = 200m','R_c = 300m')

%% play a music at the end of a code
% load chirp
% sound(y,Fs)
% load handel
% sound(y,Fs)

% lineHandles = find_system(gcs,'FindAll','On','SearchDepth',1,'Type','Line');
% Simulink.BlockDiagram.routeLine(lineHandles);