
%%%%%%%%%%%%  Script Script_samplePathDesign_circleSpiralLine.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used to generate a oval shape path ( test track ) and
% off-ramp road
% Author:       Liming
% Created Date: 2020-03-15
%
% Revisions:
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

Rc = 200; %Radius of circle arc and the end of the spiral, curvature = 1/Rc
line_length = 300; % length of the line segment, zero curvature
%1.line segment
%
lineSeg_x1 = (0:0.2:line_length)';
lineSeg_y1 = -Rc + lineSeg_x1 .* 0;
curvature_1 = zeros(length(lineSeg_y1),1);
% 2. spiral
Ls = 145;%	Length measured along the spiral curve from its initial position

[x_spiral,y_spiral,station_spiral,slope] = fnc_spiralCurve(Ls,Rc);

spiral_x1 = x_spiral' + lineSeg_x1(end);
spiral_y1 = y_spiral' + lineSeg_y1(end);
curvature_2 = slope * [0;cumsum(sqrt(diff(x_spiral').^2+diff(y_spiral').^2))];

% 3. circle
R = Rc; 
center_x2 = spiral_x1(end) - sqrt(R^2-spiral_y1(end)^2);
center_y2 = 0;

t_start = atan2((spiral_y1(end)-center_y2)/R,(spiral_x1(end)-center_x2)/R);
t= t_start:0.01:-t_start;
circleArc_x = R*cos(t);
circleArc_y = R*sin(t);

circleArc_x1 = circleArc_x'+ center_x2;
circleArc_y1 = circleArc_y'+ center_y2;
curvature_3 = (1/R).*ones(length(circleArc_x1),1);
%% 4. spiral2

spiral_x2 = flip(x_spiral') + lineSeg_x1(end);
spiral_y2 = -flip(y_spiral') + circleArc_y1(end) + max(y_spiral);

curvature_4 = slope * flip([0;cumsum(sqrt(diff(x_spiral').^2+diff(y_spiral').^2))]);


% 5. line 2 
lineSeg_x2 = (spiral_x2(end):-0.2:spiral_x2(end) - line_length)';
lineSeg_y2 = Rc + lineSeg_x2 .* 0;
curvature_5 = zeros(length(lineSeg_x2),1);
% 6. spiral3

spiral_x3 = flip(x_spiral') + lineSeg_x2(end) - max(x_spiral);
spiral_y3 = -(y_spiral') + lineSeg_y2(end) ;
curvature_6 = curvature_2;
% 7. circle2

center_x2 = spiral_x3(end) + sqrt(R^2-spiral_y3(end)^2);
center_y2 = 0;

t_start = atan2((spiral_y3(end)-center_y2)/R,(spiral_x3(end)-center_x2)/R);
t= t_start:0.01:2*pi-t_start;
circleArc_x = R*cos(t);
circleArc_y = R*sin(t);

circleArc_x2 = circleArc_x'+ center_x2;
circleArc_y2 = circleArc_y'+ center_y2;
curvature_7 = (1/R).*ones(length(circleArc_x2),1);
%% 8. spiral4

spiral_x4 = -flip(x_spiral') ;
spiral_y4 = flip(y_spiral')  + circleArc_y2(end) - max(y_spiral);

curvature_8 = curvature_4;

path.x = [lineSeg_x1(ceil(length(lineSeg_x1)/2):end);spiral_x1;circleArc_x1;spiral_x2;lineSeg_x2;spiral_x3;circleArc_x2;spiral_x4;lineSeg_x1(1:floor(length(lineSeg_x1)/2))];
path.y = [lineSeg_y1(ceil(length(lineSeg_x1)/2):end);spiral_y1;circleArc_y1;spiral_y2;lineSeg_y2;spiral_y3;circleArc_y2;spiral_y4;lineSeg_y1(1:floor(length(lineSeg_x1)/2))];
path.curvature = [curvature_1(ceil(length(lineSeg_x1)/2):end);curvature_2;curvature_3;curvature_4;curvature_5;curvature_6;curvature_7;curvature_8;curvature_1(1:floor(length(lineSeg_x1)/2))];
h_fig = figure(1);
set(h_fig,'Name','path');
x0=300;
y0=400;
width=450;
height=260;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
% plot(path.x,path.y,'k','LineWidth',2)
plot(lineSeg_x1,lineSeg_y1,'g','LineWidth',2)
plot(spiral_x1,spiral_y1,'b-.','LineWidth',2)
plot(circleArc_x1,circleArc_y1,'m:','LineWidth',2)
plot(spiral_x2,spiral_y2,'b-.','LineWidth',2)
plot(lineSeg_x2,lineSeg_y2,'g','LineWidth',2)
plot(spiral_x3,spiral_y3,'b-.','LineWidth',2)
plot(circleArc_x2,circleArc_y2,'m:','LineWidth',2)
plot(spiral_x4,spiral_y4,'b-.','LineWidth',2)

plot(path.x(1),path.y(1),'r.','MarkerSize',21)
annotation(h_fig,'textarrow',[0.579886039886041 0.517065527065528],...
    [0.321516393442623 0.235450819672131],'String',{'start'});
annotation(h_fig,'arrow',[0.591481481481483 0.502592592592594],...
    [0.802125000000002 0.803125]);
grid on
xlabel('xEast')
ylabel('yNorth')
axis equal

legend1 = legend('Line','Spiral','Arc');
set(legend1,...
    'Position',[0.28135802816081 0.440530408654682 0.152222223599752 0.165625003576279]);
% colorbar
box on



%%  remove the Duplicate data

path.station = [0;cumsum(sqrt(diff(path.x).^2+diff(path.y).^2))]; % station 

[~,unique_station_indices] = unique(path.station,'stable');
path.x  = path.x(unique_station_indices);
path.y  = path.y(unique_station_indices);
path.curvature  = path.curvature(unique_station_indices);
path.station  = path.station(unique_station_indices);

%% plot the curvature 

make_plot = 0; %parameter of fnc_parallel_curve
flag1 =1; %parameter of fnc_parallel_curve
smooth_size = 100; %parameter of fnc_parallel_curve
[~, ~, ~, ~,R_path,curvature]=fnc_parallel_curve(path.x, path.y, 1, make_plot,flag1,smooth_size);

h_fig = figure(2);
set(h_fig,'Name','path');
clf;
hold on
plot(path.station,1./R_path,'r--','LineWidth',2)

plot(path.station,path.curvature,'k.','LineWidth',2)

grid on
xlabel('station')
ylabel('curvature')
% ylim([0 0.006])
legend('calculated','Designed')
box on

% save the oval path data into *.mat file
filename = 'path.mat';
save(filename,'path')

%% off-ramp road 

h_fig = figure(889);
set(h_fig,'Name','off ramp path');
x0=300;
y0=400;
width=540;
height=200;
set(gcf,'position',[x0,y0,width,height])
clf;
hold on
% plot(path.x,path.y,'k.','LineWidth',2)
plot(lineSeg_x1,lineSeg_y1+200,'g','LineWidth',2)
plot(spiral_x1,spiral_y1+200,'b-.','LineWidth',2)

plot(circleArc_x1(1:90),circleArc_y1(1:90)+200,'c:','LineWidth',2)

plot(circleArc_x1(79),circleArc_y1(79)+200,'r.','MarkerSize',20)
grid on
xlabel('xEast')
ylabel('yNorth')
% axis equal
% colormap(jet);
% colorbar
set(gca,'XTick',[0],'YTick',[0])
xlim([0,600])
ylim([-10,180])
box on
annotation(h_fig,'textarrow',[0.231481481481481 0.27],...
    [0.453579270970576 0.287571365832236],'String',{'line'});
annotation(h_fig,'textarrow',[0.560740740740743 0.599259259259262],...
    [0.466666666666667 0.300658761528327],'String',{'spiral'});
annotation(h_fig,'textarrow',[0.826296296296299 0.79],...
    [0.388471673254283 0.493214756258235],'String',{'arc'});
annotation(h_fig,'textarrow',[0.363333333333334 0.447777777777779],...
    [0.371111111111111 0.371111111111111],...
    'Color',[0.717647058823529 0.274509803921569 1],...
    'String','Ux',...
    'LineWidth',2);

