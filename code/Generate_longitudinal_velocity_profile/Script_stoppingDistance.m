%%%%%%%%%%%%%%%%  Script_stoppingDistance.m   %%%%%%%%%%%%%%%%%%%%%
% This code is used to calculate the stopping distance for the off-ramp
% road shown as Script_samplePathDesign_circleSpiralLine.m and
% F:\Forgetting by Design Dropbox\99_Reporting and
% Publications\2021_Publications\2021_ITSC_Gao_Longitudinal Speed Preview\
% 2021_ITSC_Gao_Velocity Profile Planning with Preview_2021-04-18(with max stopping d2).docx
% 
% Author:  Liming Created Date: 2020-03-15
%
% Revisions:
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
clear all
close all

Ux0= 40; % initial speed

% Ux0= 2:0.5:50; % varying initial speed 

miu_g = 9.81*0.2; % gravity * friction coefficient 

% paramters

if length(Ux0) == 1 
    Rc_x = 10:20:440;
    Ls_y = 5:20:440;
else
    Rc_x = 200;
    Ls_y = 145;
end

[Rc,Ls] = meshgrid(Rc_x,Ls_y);

%% calculate the stopping distance backward
%% start from: arc segment
Dsa = pi*Rc./4; %maximum stoooing distance at circular arc
Uxa0 = sqrt(miu_g.*Rc); % maximum initial speed at arc 

%% then: spiral segment 

slope  = 1./(Rc.*Ls);

N = 1000; % the numcer of integral steps
station_span = zeros(length(Ls),N );
for i = 1:length(Ls)
    station_span(i,:) = linspace(0,Ls(i),N );
end

% Ux_curceLimit = sqrt(miu_g ./(slope .* (Ls-station_span)) );

% % ode method
% f = @(s,Ux)(1/Ux)* sqrt(miu_g^2 - (slope*(Ls-s)*Ux^2)^2);
% [s,Ux_ode] = ode45(f, station_span,Ux_0);
% Ux_ode_real = real(Ux_ode);

%  forward integration
station_backInt = flip(station_span,2);
delta_s = -diff(station_backInt,1,2);

Ux_backInt = zeros([size(station_span),length(Rc_x)]);

for i_Uxa0 = 1:length(Rc_x)
    
    Ux_backInt(:,1,i_Uxa0) = Uxa0(i_Uxa0);
    
    for i = 1:length(station_backInt)-1
        Ux_backInt(:,i+1,i_Uxa0) = Ux_backInt(:,i,i_Uxa0) + (delta_s(:,i)./Ux_backInt(:,i,i_Uxa0)).* sqrt(miu_g.^2 - (slope(:,i_Uxa0).*(Ls(:,i_Uxa0)-station_span(:,i)).*Ux_backInt(:,i,i_Uxa0).^2).^2);
    end
    
end

%% plot
h_fig = figure(210);
set(h_fig,'Name','path');
clf;
hold on
% plot(station_span,Ux_curceLimit,'b','LineWidth',2)
% plot(station_span,Ux_ode_real,'g','LineWidth',2)
for i_Uxa0 = 1:length(Rc_x)
    plot(station_span,Ux_backInt(:,:,i_Uxa0),'m.','LineWidth',2)
end
grid on
xlabel('station')
ylabel('velocity (m/s)')
ylim([0 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on

%% finally: line segment
% Ux_fl = max(Ux_ode_real)
Ux_fl = zeros(size(slope));
for i_Uxa0 = 1:length(Rc_x)
    aa= squeeze(Ux_backInt(:,:,i_Uxa0));
    Ux_fl(:,i_Uxa0) = real(max(aa,[],2)); %
end

D_sl = (Ux0.^2 - Ux_fl.^2)./(2*miu_g);


%% add the distance at each segment 
if length(Ux0) == 1
    D_total  = Dsa + Ls + D_sl;
else
    D_total = zeros(size(Ux0));
    for i = 1: length(D_sl)
        if Ux0(i) < min(Ux_backInt)
            D_total(i) = (Rc_x/2)* atan(sqrt(Ux0(i)^4/(miu_g^2*Rc_x^2 -Ux0(i)^4 )))
        elseif Ux0(i) > max(Ux_backInt)
            D_total(i)  = Dsa + Ls + D_sl(i);
        else
            [Idx,D] = knnsearch(Ux_backInt',Ux0(i));
            Ls_empty = station_span(length(station_span)+1 - Idx);
            D_total(i)  = Dsa + Ls -Ls_empty ;
        end
            
    end
end

%% plot the results
if length(Ux0) == 1
    h_fig = figure(2100);
    set(h_fig,'Name','path');
    x0=300;
    y0=400;
    width=540;
    height=420;
    set(gcf,'position',[x0,y0,width,height])
    clf;
    hold on
    
    s = surf(Rc,Ls,D_total,'FaceAlpha',0.6);
    contour(Rc,Ls,D_total,'ShowText','on')
    s.EdgeColor = 'none';
    % contour3(Rc,Ls,D_total)
    view(-43.5,19)
    colormap jet
    %plot(Rc,D_total,'b.','LineWidth',2)
    % plot(station_span,Ux_ode_real,'g','LineWidth',2)
    % plot(station_span,Ux_backInt,'m--','LineWidth',2)
    ylim([0 430])
    xlim([0 430])
    grid on
    xlabel('Rc (m)')
    ylabel('Ls (m)')
    zlabel('Preview Distance(m)')
    % ylim([0 60])
    % axis equal
    % legend('Curve limit speed','vehicle dynamics limit speed')
    box on
     
    h_fig = figure(21200);
    set(h_fig,'Name','path');
    x0=300;
    y0=400;
    width=540;
    height=420;
    % set(gcf,'position',[x0,y0,width,height])
    clf;
    hold on
    
    s = surf(1./Rc,slope,D_total,'FaceAlpha',1);
    view(-39,36)
    colormap jet
    %plot(Rc,D_total,'b.','LineWidth',2)
    % plot(station_span,Ux_ode_real,'g','LineWidth',2)
    % plot(station_span,Ux_backInt,'m--','LineWidth',2)
    
    grid on
    xlabel('Curvature (1/m)')
    ylabel('Ls (m)')
    zlabel('Stop Distance(m)')
    % ylim([0 60])
    % axis equal
    % legend('Curve limit speed','vehicle dynamics limit speed')
    box on
else
    h_fig = figure(2100);
    set(h_fig,'Name','path');
    x0=300;
    y0=400;
    width=540;
    height=300;
    set(gcf,'position',[x0,y0,width,height])
    clf;
    hold on
    plot(Ux0,D_total,'b-','LineWidth',2)
    plot(min(Ux_backInt)*ones(1,1000),[1:1:1000],'k-.','LineWidth',1)
    plot(max(Ux_backInt)*ones(1,1000),[1:1:1000],'k-.','LineWidth',1)
%     ylim([0 420])
%     xlim([0 420])
    grid on
    xlabel('U_x_0 (m/s)')
    ylabel('Preview Distance (m)')

    % axis equal
    % legend('Curve limit speed','vehicle dynamics limit speed')
    box on
    
end