%%%%%%%%%%%%  Script Script_speedProfile_samplePath.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Purpose: Generate the limit speed profile given a path 

% Author:       Liming
% Created Date: 2021-04-01
%
% Reference: 

%
% To do list:
%     need a better version of piecewise segment algorithm 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: Prepare the workspace and flag
% Clear workspace, figures, and console
clear all;  %#ok<CLALL> % Clears workspace to remove any old variables
close all;  % Close all open figures
clc;        % Clear console space

flag.friction = true; % true: use changing friction 
flag.grade = false; % true: use non-zero grade 
flag.curvature = true; %true: use non-zero curvature 

%%Initialize parameters that may change
lamuda = 0.95;
speed_limit = 50; %m/s
default_friction = 0.1; % used for constant friction assumption

g = 9.81; %m/s^2

% for path_desired, lamuda = 0.4, speed_limit = 25
%% step 2: load a path including position and curvature, and then assign hypothetical friction and grade
load('path.mat','path');

path_station = max(path.station);

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

% assign friction distribution to the path  
friction_coefficient  = zeros(length(path.station),1);
% high friction
index_fhigh = ((path.station < 0.2*path_station) | (path.station > 0.85*path_station));
friction_coefficient(index_fhigh)  = 0.9; 
index_f1 = ( path.station <= 0.3750*path_station & path.station >= 0.200*path_station);
friction_coefficient(index_f1)  = 0.3; 
index_f2 = ( path.station <= 0.700*path_station & path.station >= 0.550*path_station);
friction_coefficient(index_f2)  = 0.15; 
index_f3 = ( path.station <= 0.850*path_station & path.station >= 0.700*path_station);
friction_coefficient(index_f3)  = 0.35; 
index_f4 = ( path.station < 0.550 * path_station & path.station > 0.3750*path_station);
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

% assign grade distribution 
grade  = zeros(length(path.station),1);
if flag.grade
    index_grade1 = ((path.station < 0.15*path_station) & (path.station > 0.05*path_station));
    grade(index_grade1)  = 6;  % grade unit: percentage
    
    index_grade2 = ( path.station < 0.50 * path_station & path.station > 0.40*path_station);
    grade(index_grade2)  = -6;
    
end

path.grade = atan(grade./100); % percentage to radii

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

%% step2.5 cut the path into segments by curvature, friction, and grade.(need a better version)
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
%% step 2: plan velocity profile 

% 1. find the curve speed limit of the road given the friction and geomtry
% 

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
plot(path.station(index_bound),U_max(index_bound),'k.','MarkerSize',20)

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

ay = v.^2.*path.curvature';
ax = [sqrt((path.friction(1)*g)^2 - ay(1).^2) 0.5.*diff(v.^2)./diff(path.station')];

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
annotation(h_fig,'textbox',...
    [0.811370370370371 0.157333333333333 0.0893703703703704 0.0906666666666668],...
    'String','(a)',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.19803703703704 0.616000000000004 0.0893703703703703 0.0906666666666668],...
    'String','A',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.236555555555559 0.368000000000003 0.0893703703703703 0.0906666666666665],...
    'String','B',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.391370370370374 0.820000000000004 0.0893703703703703 0.0906666666666661],...
    'String','C''',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.44618518518519 0.574666666666672 0.089370370370371 0.0906666666666661],...
    'String','C',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.573592592592598 0.273333333333335 0.0893703703703704 0.0906666666666668],...
    'String','D',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
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
filename = 'path_desired.mat'; % 
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
plot(200,35,'m.','MarkerSize',20)

grid on
xlabel('Station (m)')
ylabel('Ux (m/s)')
xlim([0 max(path.station)])
ylim([0 speed_limit+5])
legend('Curve limit speed','Speed profile for actual ¦Ì','Speed profile for safe')

% legend('Curve limit speed','Speed profile for actual ¦Ì','Speed profile for confident low bound ¦Ì = 0.2')
box on
%% calculate the preview distance 
% v = min(real(Ux_forwardInt), real(Ux_backwardInt));
% difference_max_profile = U_max - v';  % thid may not be the correct preview distance 
Ux_backwardInt = real(Ux_backwardInt);
difference_max_profile = U_max - Ux_backwardInt';
index_plan = find(difference_max_profile>0.1);
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
% plot(900.4,50,'k.','MarkerSize',25) % C'
plot(900.4,445,'r.','MarkerSize',25) % C
plot(1345.5,0,'r.','MarkerSize',25 ) % D
plot(419,0,'r.','MarkerSize',25)% B
plot(315,110,'r.','MarkerSize',25)% A

annotation(h_fig,'textbox',...
    [0.252851851851855 0.340800000000004 0.0893703703703702 0.0906666666666668],...
    'String','A',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.278037037037041 0.185600000000003 0.0893703703703703 0.0906666666666665],...
    'String','B',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
    [0.474333333333339 0.828000000000001 0.0893703703703709 0.0981333333333381],...
    'String','C',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(h_fig,'textbox',...
     [0.622481481481487 0.186933333333335 0.0893703703703704 0.0906666666666668],...
    'String','D',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');

box on
%%
d_preview_01_50 = d_preview;
% save('d_preview_01_50.mat','d_preview_01_50')
h_fig = figure(455078);
set(h_fig,'Name','diff');
x0=300;
y0=400;
width=540;
height=250;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
load('d_preview_08_45.mat')
load('d_preview_05_45.mat')
load('d_preview_02_45.mat')
load('d_preview_02_50.mat')
load('d_preview_01_50.mat')
% plot(path.station,d_preview,'m','LineWidth',2)
plot(path.station,d_preview_02_45,'g-','LineWidth',1.5)
plot(path.station,d_preview_05_45,'b-.','LineWidth',2)
plot(path.station,d_preview_08_45,'r:','LineWidth',2)
grid on
xlabel('Station (m)')
ylabel('Preview Distance (m)')
xlim([0 max(path.station)])
% ylim([0 speed_limit+5])

legend('¦Ì_c_o_n_f = 0.2','¦Ì_c_o_n_f = 0.5','¦Ì_c_o_n_f = 0.8')
box on

%
index_dangerous_true = find(d_preview>0) ;
index_dangerous_01 = find(d_preview_01_50>0) ;
h_fig = figure(478978);
set(h_fig,'Name','diff');
x0=300;
y0=400;
width=540;
height=300;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
plot(path.station,path.curvature,'b','LineWidth',2)
plot(path.station(index_dangerous_01),path.curvature(index_dangerous_01),'go','LineWidth',2)
plot(path.station(index_dangerous_true),path.curvature(index_dangerous_true),'r.','LineWidth',2)
grid on
xlabel('station (m)')
ylabel('curvature (m/s)')
% xlim([0 max(path.station)])
% ylim([0 speed_limit+5])
legend('all curvature','Dangerous station for confident low bound ¦Ì = 0.1','Dangerous station for actual ¦Ì > 0.15')
box on
%

%% acceleration
h_fig = figure(455022);
set(h_fig,'Name','acceleariton');
x0=300;
y0=400;
width=540;
height=300;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
plot(path.station,ax,'b-.','LineWidth',1.5)
plot(path.station,ay,'g--','LineWidth',1.5)
plot(path.station,sqrt(ax.^2 + ay.^2),'r:','LineWidth',2)
plot(path.station,a_max,'k-','LineWidth',1.5)
grid on
xlabel('Station(m)')
ylabel('Acceleration (m/s^2)')
xlim([0 max(path.station)])
% ylim([10 speed_limit])
legend1 = legend('$$a_x$$','$a_y$','$$a_{total}$$','$$a_{max}$$');

% h = legend(['$$\sqrt{blah}$$'])
% set(h,'Interpreter','latex','fontsize',24) 

% 'sqrt()'
set(legend1,...
    'Interpreter','latex',...
    'FontSize',11,...
    'Position',[0.348148148148149 0.143333333333334 0.290567938919719 0.163611110983955],...
    'NumColumns',2);

annotation(h_fig,'textbox',...
    [0.811370370370371 0.157333333333333 0.0893703703703704 0.0906666666666668],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');

box on


h_fig = figure(49822);
set(h_fig,'Name','acceleariton');
x0=300;
y0=400;
width=400;
height=400;
% set(gcf,'position',[x0,y0,width,height])
clf;
hold on
plot(ay,ax,'b.','LineWidth',2)
xlim([-10 10])
ylim([-10 10])
grid on
xlabel('ay (m/s^2)')
ylabel('ax (m/s^2)')
axis equal
box on