function[] = OverheadPlot(figHandle, t, posE, posN, psi, delta, plotTimes, vehicle)
% Function that plots an overhead view of the car trajectory for a given dataset
% Modified for ME227 by Craig Beal
% Date: 5.26.2010
%
% OverheadPlot is a function that allows you to plot an overhead view of the car
% in global coordinates, including the steering angles on each wheel.  Call
% OverheadPlot like:
%
% OverheadPlot(figHandle, t, posE, posN, psi, delta, plotTimes, vehicle);
%
% figHandle: is a handle to a figure; get this like figHandle = figure; or
%       figHandle = figure(1);
% t: your simulation time vector
% posE and posN: the vectors of your global east and north positions
% psi: your global heading vector
% delta: a matrix of steering angles with 4 angles at each time step
% plotTimes: a vector you define to tell the plot code which time steps to
%       draw the vehicle body
% vehicle: your structure of vehicle parameters (a, b, d, and Re are needed)

% Check to make sure that the data vectors are the same length
N = length(t);
if length(posE) ~= N || length(posN) ~= N || length(psi) ~= N
    error('One or more data vectors differ in length');
end

% Set the total number of positions of the car to plot
Ncars = length(plotTimes);

% Check the plotTimes vector to make sure all the times are in the data
for i = 1:length(Ncars)
    if plotTimes(i) < t(1) || plotTimes(i) > t(end)
        error('One or more entries in plotTimes outside your time vector data');
    end
end

% Vehicle Parameters
a = vehicle.a;
b = vehicle.b;                   
t_f = vehicle.d;
t_r = vehicle.d;
r_w = vehicle.Re;

plotindex = zeros(Ncars,1);
for i = 1:Ncars
    % get index corresponding to plotTimes
    delT = t(2) - t(1);
    plotindex(i) = floor(plotTimes(i)/delT+1);
end

figure(figHandle);
hold on;
axis equal
grid on;

% Check to see that delta is of the right size
[m,n] = size(delta);
if max(m,n) ~= N
    error('Steering data matrix does not match time reference');
end
if ~(m==4 || n==4)
    error('Wrong number of steering angles per time step');
else
    if m==4
        delta = delta';
    end
end

for i = 1:Ncars
    PosEplot = posE(plotindex(i));
    PosNplot = posN(plotindex(i));
    Heading = psi(plotindex(i));
    Delta = delta(plotindex(i),:);
    DrawVehicle(PosEplot,PosNplot,Heading,Delta,a,b,t_f,t_r,r_w);
end
% plot trajectory of vehicle CG
% figure(8); hold on;
plot(posE(plotindex(1):plotindex(end)),posN(plotindex(1):plotindex(end)),'k')

axes_handle = gca;
set(axes_handle, 'FontName', 'Times New Roman', 'FontSize', 14);
title('Vehicle Trajectory Top View','FontSize', 14);
xlabel('East Position (m)', 'FontSize', 14)
ylabel('North Position (m)', 'FontSize', 14)


% A helper function for the OverheadPlot script
function DrawVehicle(PosE, PosN, Heading, Delta, a, b, t_f, t_r, r_w)

lf = 1; rf = 2; lr = 3; rr = 4;

%Draw line from CG to front axle, CG to rear axle - the "body"

Heading = Heading - pi/2;

FrontAxle_Center_x = PosE - a*sin(Heading);
FrontAxle_Center_y = PosN + a*cos(Heading);

RearAxle_Center_x = PosE + b*sin(Heading);
RearAxle_Center_y = PosN - b*cos(Heading);

FrontBody = [PosE FrontAxle_Center_x; PosN FrontAxle_Center_y];
RearBody=  [PosE RearAxle_Center_x; PosN RearAxle_Center_y];

plot(PosE, PosN, 'ko', 'MarkerSize', 4)
plot1=plot(FrontBody(1,:), FrontBody(2,:), 'LineWidth', 2.5);
set(plot1(1),'Color',[.8 .8 .8]);
plot1=plot(RearBody(1,:), RearBody(2,:), 'LineWidth', 2.5);
set(plot1(1),'Color',[.8 .8 .8])

%Now draw axles
FrontAxle_Right_x = FrontAxle_Center_x + (t_f/2)*cos(Heading);
FrontAxle_Right_y = FrontAxle_Center_y + (t_f/2)*sin(Heading);

FrontAxle_Left_x = FrontAxle_Center_x - (t_f/2)*cos(Heading);
FrontAxle_Left_y = FrontAxle_Center_y - (t_f/2)*sin(Heading);

RearAxle_Right_x = RearAxle_Center_x + (t_r/2)*cos(Heading);
RearAxle_Right_y = RearAxle_Center_y + (t_r/2)*sin(Heading);

RearAxle_Left_x = RearAxle_Center_x - (t_r/2)*cos(Heading);
RearAxle_Left_y = RearAxle_Center_y - (t_r/2)*sin(Heading);

FrontAxle = [FrontAxle_Left_x FrontAxle_Right_x; FrontAxle_Left_y FrontAxle_Right_y];
RearAxle = [RearAxle_Left_x RearAxle_Right_x; RearAxle_Left_y RearAxle_Right_y];

plot1=plot(FrontAxle(1,:), FrontAxle(2,:), 'LineWidth', 2.5);
set(plot1(1),'Color',[.8 .8 .8])
plot1=plot(RearAxle(1,:), RearAxle(2,:), 'LineWidth', 2.5);
set(plot1(1),'Color',[.8 .8 .8])

%Now draw wheels

RightFrontTire_Front_x = FrontAxle_Right_x - r_w*sin(Heading+Delta(rf));
RightFrontTire_Front_y = FrontAxle_Right_y + r_w*cos(Heading+Delta(rf));

RightFrontTire_Back_x = FrontAxle_Right_x + r_w*sin(Heading+Delta(rf));
RightFrontTire_Back_y = FrontAxle_Right_y - r_w*cos(Heading+Delta(rf));

RightRearTire_Front_x = RearAxle_Right_x - r_w*sin(Heading+Delta(rr));
RightRearTire_Front_y = RearAxle_Right_y + r_w*cos(Heading+Delta(rr));

RightRearTire_Back_x = RearAxle_Right_x + r_w*sin(Heading+Delta(rr));
RightRearTire_Back_y = RearAxle_Right_y - r_w*cos(Heading+Delta(rr));

LeftFrontTire_Front_x = FrontAxle_Left_x - r_w*sin(Heading+Delta(lf));
LeftFrontTire_Front_y = FrontAxle_Left_y + r_w*cos(Heading+Delta(lf));

LeftFrontTire_Back_x = FrontAxle_Left_x + r_w*sin(Heading+Delta(lf));
LeftFrontTire_Back_y = FrontAxle_Left_y - r_w*cos(Heading+Delta(lf));

LeftRearTire_Front_x = RearAxle_Left_x - r_w*sin(Heading+Delta(lr));
LeftRearTire_Front_y = RearAxle_Left_y + r_w*cos(Heading+Delta(lr));

LeftRearTire_Back_x = RearAxle_Left_x + r_w*sin(Heading+Delta(lr));
LeftRearTire_Back_y = RearAxle_Left_y - r_w*cos(Heading+Delta(lr));


RightFrontTire = [RightFrontTire_Front_x RightFrontTire_Back_x;... 
    RightFrontTire_Front_y RightFrontTire_Back_y];
RightRearTire = [RightRearTire_Front_x RightRearTire_Back_x;... 
    RightRearTire_Front_y RightRearTire_Back_y];
LeftFrontTire = [LeftFrontTire_Front_x LeftFrontTire_Back_x;... 
    LeftFrontTire_Front_y LeftFrontTire_Back_y];
LeftRearTire = [LeftRearTire_Front_x LeftRearTire_Back_x;... 
    LeftRearTire_Front_y LeftRearTire_Back_y];

plot(RightFrontTire(1,:), RightFrontTire(2,:), 'r', 'LineWidth', 3)
plot(RightRearTire(1,:), RightRearTire(2,:), 'k', 'LineWidth', 3)
plot(LeftFrontTire(1,:), LeftFrontTire(2,:), 'r', 'LineWidth', 3)
plot(LeftRearTire(1,:), LeftRearTire(2,:), 'k', 'LineWidth', 3)
