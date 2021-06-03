%%%%%%%%%%%%%%%%  Script_spiralPathDesign.m   %%%%%%%%%%%%%%%%%%%%%
% This code is used to:
% Section 1: test code to generate the spiral curve. The conclusion is used at fnc_spiralCurve 
% Section 2: given the paramers of a spiral, calculate the acceleration behavior
% section 3: given the paramers of a spiral, calculate the stopping distance 
% Author:  Liming Created Date: 2020-03-15
%
% Revisions:
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%% ========== section1: plot the x and y coordinates given the parameters of a spiral =======

%% Normallized Euler curve %https://en.wikipedia.org/wiki/Euler_spiral
Ls = 145;%	Length measured along the spiral curve from its initial position
delta_s = 0.1;
% Rc = 1/(2*Ls); %Radius of circular curve at the end of the spiral
t = 0:delta_s:45;
len = length(t);
x_spiral_normal = zeros(1, len);
y_spiral_normal = zeros(1, len);
fun_c = @(s) cos(s.^2);
fun_s = @(s) sin(s.^2);

for n = 2:len
    x_spiral_normal(n) = x_spiral_normal(n-1)+ integral(fun_c, t(n-1), t(n));
    y_spiral_normal(n) = y_spiral_normal(n-1)+ integral(fun_s, t(n-1), t(n));
    
end

% q = integral(@(s)fun_c(s),0,2.5);
make_plot = 0; %parameter of fnc_parallel_curve
flag1 =0; %parameter of fnc_parallel_curve
smooth_size = 0; %parameter of fnc_parallel_curve
[~, ~, ~, ~,R_spiral_normal,UnitNormalV]=fnc_parallel_curve(x_spiral_normal, y_spiral_normal, 1, make_plot,flag1,smooth_size); % calculate the curvature
station_spiral_normal = [0 cumsum(sqrt(diff(x_spiral_normal).^2+diff(y_spiral_normal).^2))]; % calculate the station
% curvature  = 2*station 

%% denormalize Euler curve
Rc = 200; %Radius of circular curve at the end of the spiral

%LS = 45;
factor_a= 1/sqrt(2*Rc*Ls); % a constant scaling up (denormalize) factor 
slope = 2*factor_a^2
% curvatre = 2*factor_a^2 * Ls
s_spiral = 0:delta_s*factor_a:Ls*factor_a; % station

x_spiral = zeros(1, length(s_spiral));
y_spiral = zeros(1, length(s_spiral));
for n = 2:length(s_spiral)
    x_spiral(n) = x_spiral(n-1)+ integral(fun_c, s_spiral(n-1), s_spiral(n))/factor_a;
    y_spiral(n) = y_spiral(n-1)+ integral(fun_s, s_spiral(n-1), s_spiral(n))/factor_a;
    
end

% q = integral(@(s)fun_c(s),0,2.5);
[~, ~, ~, ~,R_spiral,UnitNormalV]=fnc_parallel_curve(x_spiral, y_spiral, 1, make_plot,flag1,smooth_size);
station_spiral = [0 cumsum(sqrt(diff(x_spiral).^2+diff(y_spiral).^2))];
curvature  = 1./R_spiral;

%% cubic fit 
x = 0:0.1:Ls;
% y = sqrt(100-x.^2);
y = x.^3/(6*Rc*Ls);
station_span = [0 cumsum(sqrt(diff(x).^2+diff(y).^2))];

[~, ~, ~, ~,R,UnitNormalV]=fnc_parallel_curve(x, y, 1, make_plot,flag1,smooth_size);
    

%% plot
h_fig = figure(1);
set(h_fig,'Name','path');
clf;
hold on
% plot(x_spiral_normal,y_spiral_normal,'b','LineWidth',2)
plot(x_spiral,y_spiral,'k','LineWidth',2)
% plot(x,y,'r--','LineWidth',2)
% scatter3(path.center_x,path.center_y,path.friction,10, path.friction,'.');
grid on
xlabel('xEast')
ylabel('yNorth')
% axis equal
% colormap(jet);
% colorbar
box on
% ylim([-5 5])

h_fig = figure(2);
set(h_fig,'Name','path');
clf;
hold on
% plot(station,1./R,'r--','LineWidth',2)
% plot(station_spiral_normal,1./R_spiral_normal,'b','LineWidth',2)
plot(station_spiral(1:end-2),curvature(1:end-2),'k','LineWidth',2)
% scatter3(path.center_x,path.center_y,path.friction,10, path.friction,'.');
grid on
xlabel('station')
ylabel('curvature')
% axis equal
box on

%% try to the analytical solution for speed equation, but failed 
syms y(t) a b  x c
ode = (diff(y,t)+y)^2 == 1;
cond = y(0) == 0;
ySol(t) = dsolve(ode,cond)


ode2 = (diff(y,t))^2 == a^2- y^4*b*t^2;
cond2 = y(0) == 10;
ySol2(t) = dsolve(ode2,cond2)

ode3 = diff(y,t) == (1/y)*sqrt(a - b*y^4);
cond3 = y(0) == 10;
ySol3(t) = dsolve(ode3,cond3)

ode4 = diff(y,t) == sqrt(a - b*y^4);
cond4 = y(0) == 10;
ySol4(t) = dsolve(ode4,cond4)

y1 = x/(c+sqrt(a^2-b^2*x^4));
fy1=int(y1)
symvar(y1,1)

%% ========== section 2: given the paramers of a spiral, calculate the acceleration behavior =======
%% 
Ls = 150 % spiral length
Rc = 200 % radius at the end of spiral 
slope  = 1/(Rc*Ls) % curvature = slope * station
theta= 0*5.71*pi/180; % grade

g= 9.81; % gravity acceleration
miu_g = 0.8 * g; % ¦Ì*g, ¦Ì is friction coefficient

Ux_0 = sqrt(miu_g /(slope * Ls));  % maximum allowable inital speed 
% Ux_0 = 10; % other initial speed 
Ux_max = 60; % the road traffic speed limit 


station_span = 0:0.01:Ls;
curvature_span = station_span * slope;
curvature_span_back = flip(curvature_span);

Ux_curceLimit = sqrt(miu_g ./(slope .* (Ls-station_span)) ); % curve limit speed

%% method 1: numerical soultion

% == numerical solution 1: ode 45  (it is more accurate )
f = @(s,Ux)(1/Ux)* (sqrt((miu_g*cos(theta))^2 - (slope*(Ls-s)*Ux^2)^2) +  g*sin(theta)); % acceleration 

[s,Ux_frictionCircle] = ode45(f, station_span,Ux_0);
Ux_frictionCircle_real = real(Ux_frictionCircle);

% == numerical solution 2: backward numerical integration (thisi based on the Euler method)

delta_s = diff(station_span);

Ux_backInt = zeros(1,length(station_span));
Ux_backInt(1) = Ux_0;
for i = 1:length(delta_s)
    Ux_backInt(i+1) = Ux_backInt(i) + (delta_s(i)/Ux_backInt(i))* (sqrt((cos(theta)*miu_g)^2 - (slope*(Ls-station_span(i))*Ux_backInt(i)^2)^2)+ g*sin(theta));
end

h_fig = figure(210);
set(h_fig,'Name','Ux ');
clf;
hold on
plot(station_span,Ux_curceLimit,'b','LineWidth',2)
plot(station_span,Ux_frictionCircle_real,'g','LineWidth',2)
plot(station_span,Ux_backInt,'b--','LineWidth',2)

Ux_backIntddd = flip(Ux_backInt);
grid on
xlabel('station')
ylabel('velocity (m/s)')
ylim([0 Ux_max])
% axis equal
legend('Curve limit speed',' friction limit speed ode45','friction limit speed euler','Location','best')
box on

%% method 2 piecewise circle arc approximation (get a conservative results)

% find the station where curve speed limit reaches the maximum speed limit
% curvature_Ux_max = miu_g / Ux_max^2;
% Idx_cur = knnsearch(curvature_span',curvature_Ux_max); 

% fit spiral using arcs

nbs_arcs = 10; % number of arcs for approximation
% nbs_arcs = ceil(Ls/(Rc*pi/4))

% even piecewise
index_step  = floor(length(station_span)/nbs_arcs); % step size to cut the spiral into nbs_arcs pieces

curvature_arcs = [curvature_span_back(1:index_step:index_step*(nbs_arcs-1)+1) curvature_span_back(index_step*(nbs_arcs-1)+1)]; % start curavtures for each pieces, the last one is the curvature of residual segment 
curvature_arcs_span_back = curvature_span_back;
for i= 1:length(curvature_arcs)
    
    index_end = min(index_step*i,length(curvature_arcs_span_back));
    curvature_arcs_span_back(index_step*(i-1)+1:index_end) = curvature_arcs(i);
    
end

% curvature_arcs_span_back(index_step*length(curvature_arcs)+1:end) = curvature_arcs(end);

Ux_arc_curceLimit = sqrt(miu_g ./curvature_arcs_span_back );
Ux_arc_curceLimit(Ux_arc_curceLimit>Ux_max) = Ux_max;

% accelerate piecewise using closed form arc solution
Ux_arc_approximation = zeros(size(Ux_arc_curceLimit));
u_x_initial = Ux_0;
for i= 1:length(curvature_arcs)
    index_end = min(index_step*i,length(curvature_arcs_span_back));
    c1 = (1/(2*curvature_arcs(i)))* atan(sqrt(u_x_initial^4/((miu_g/curvature_arcs(i))^2 -u_x_initial^4 )));
    station_c = station_span(index_step*(i-1)+1:index_end) - station_span(index_step*(i-1)+1);
    u_x_acce = sqrt(miu_g/curvature_arcs(i)).* sqrt(tan(2*(c1+station_c).* curvature_arcs(i)))./((tan(2.*(c1+station_c).*curvature_arcs(i)).^2+1).^(1/4)); % the period of this equation is pi*R/2
    u_x_acce(real(u_x_acce)<=0) = max(real(u_x_acce)); % keep speed
    
    Ux_arc_approximation(index_step*(i-1)+1:index_end) = u_x_acce;
    u_x_initial = u_x_acce(end);
    clear c1 u_x_acce
end

max_error = max(Ux_frictionCircle_real' - Ux_arc_approximation)


% arithmetic piecewise
seg_initial = 500;
increasement_d = floor(2*(index_step*nbs_arcs - nbs_arcs*seg_initial)/(nbs_arcs*(nbs_arcs-1)));

seg_interval = seg_initial: increasement_d: nbs_arcs*(increasement_d-1)+seg_initial;

index_curvature = [0 cumsum(seg_interval) length(curvature_span_back)];

curvature_arcs_arith = [curvature_span_back(index_curvature(1:end-2)+1) curvature_span_back(index_curvature(end-2)+1)]; % start curavtures for each pieces, the last one is the curvature of residual segment 
curvature_arcs_span_back_arith = curvature_span_back;
for i= 1:length(curvature_arcs_arith)
    
    index_end = min(index_curvature(i+1),length(curvature_arcs_span_back_arith));
    index_start =index_curvature(i)+1;
    curvature_arcs_span_back_arith(index_start:index_end) = curvature_arcs_arith(i);
end

Ux_arc_approximation_arith = zeros(size(Ux_arc_curceLimit));
u_x_initial = Ux_0;
for i= 1:length(curvature_arcs_arith)
    index_start =index_curvature(i)+1;
    index_end = min(index_curvature(i+1),length(curvature_arcs_span_back_arith));
    
    c1 = (1/(2*curvature_arcs_arith(i)))* atan(sqrt(u_x_initial^4/((miu_g/curvature_arcs_arith(i))^2 -u_x_initial^4 )));
    station_c = station_span(index_start:index_end) - station_span(index_start);
    u_x_acce = sqrt(miu_g/curvature_arcs_arith(i)).* sqrt(tan(2*(c1+station_c).* curvature_arcs_arith(i)))./((tan(2.*(c1+station_c).*curvature_arcs_arith(i)).^2+1).^(1/4)); % the period of this equation is pi*R/2
    u_x_acce(real(u_x_acce)<=0) = max(real(u_x_acce)); % keep speed
    
    Ux_arc_approximation_arith(index_start:index_end) = u_x_acce;
    u_x_initial = u_x_acce(end);
    clear c1 u_x_acce
end

max_error = max(Ux_frictionCircle_real' - Ux_arc_approximation_arith)

h_fig = figure(21145);
set(h_fig,'Name','Curvature ');
clf;
hold on
plot(station_span,flip(curvature_arcs_span_back),'b','LineWidth',2)
plot(station_span,curvature_span,'g','LineWidth',2)
plot(station_span,flip(curvature_arcs_span_back_arith),'r','LineWidth',2)
grid on
xlabel('station')
ylabel(',curvature (1/m)')
% ylim([0 Ux_max+10])
% axis equal
legend('curvature arcs','curvature spiral','Location','best')
box on


h_fig = figure(2170);
set(h_fig,'Name','Ux ');
clf;
hold on
plot(station_span,Ux_curceLimit,'b','LineWidth',2)
plot(station_span,Ux_arc_curceLimit,'m','LineWidth',2)
plot(station_span,Ux_frictionCircle_real,'g','LineWidth',2)
plot(station_span,Ux_arc_approximation,'b--','LineWidth',2)
plot(station_span,Ux_arc_approximation_arith,'r--','LineWidth',2)
grid on
xlabel('station')
ylabel('velocity (m/s)')
ylim([0 Ux_max+10])
% axis equal
legend('Actual Curve limit speed','Approximate Curve limit speed','friction limit speed ode45','friction limit arcs approximation','Location','best')
box on

h_fig = figure(7470);
set(h_fig,'Name','Ux ');
clf;
hold on
plot(station_span,Ux_frictionCircle_real,'g','LineWidth',2)
plot(station_span,Ux_arc_approximation,'b--','LineWidth',2)
plot(station_span,Ux_arc_approximation_arith,'r--','LineWidth',2)
grid on
xlabel('station')
ylabel('velocity (m/s)')
% ylim([0 Ux_max])
% axis equal
legend('friction limit speed ode45','friction limit even arcs approximation','friction limit arithmetic arcs approximation','Location','best')
box on


%% method 3 Taylor's series approximation solution ( polynominal soultion)


syms y(x) a b x c d s gs
assume(a >= 0)
dece = false; % deceleration or acceleration
if dece
    y_p = -(1/y)*sqrt(a^2-(b*x*y^2)^2);
    
else
    y_p = (1/y)*(sqrt(a^2-(b*(d-x)*y^2)^2)+gs);
    
end

y_p_zero = subs(y_p,[x y],[0 c]);

% G = functionalDerivative(f,y)
y_pp = diff(y_p,x);

y_ppp = diff(y_pp,x); % 3rd order derivative 

y_pppp = diff(y_ppp,x);

y_5p = diff(y_pppp,x);

y_6p = diff(y_5p,x);
% argnames(y_pp)

% x = 0

y_pp_imd = subs(y_pp,[x diff(y(x), x)],[0 y_p_zero]);
y_pp_zero = subs(y_pp_imd,y(0),c);
% simplify(y_pp_zero);

y_ppp_imd = subs(y_ppp,[x diff(y(x), x) diff(y(x), x, x)],[0 y_p_zero y_pp_zero]);
y_ppp_zero = simplify(subs(y_ppp_imd,y(0),c));
% simplify(y_ppp_zero);


y_pppp_imd = subs(y_pppp,[x diff(y(x), x) diff(y(x), x, x) diff(y(x), x, x,x)],[0 y_p_zero y_pp_zero y_ppp_zero]);
y_pppp_zero = simplify(subs(y_pppp_imd,y(0),c));
% simplify(y_pppp_zero);


y_5p_imd = subs(y_5p,[x diff(y(x), x) diff(y(x), x, x) diff(y(x), x, x,x) diff(y(x), x, x,x,x)],[0 y_p_zero y_pp_zero y_ppp_zero y_pppp_zero]);
y_5p_zero = subs(y_5p_imd,y(0),c);
simplify(y_5p_zero);

y_6p_imd = subs(y_6p,[x diff(y(x), x) diff(y(x), x, x) diff(y(x), x, x,x) diff(y(x), x, x,x,x) diff(y(x), x, x,x,x,x)],[0 y_p_zero y_pp_zero y_ppp_zero y_pppp_zero y_5p_zero   ]);
y_6p_zero = subs(y_6p_imd,y(0),c);
simplify(y_6p_zero);

coefficient = [c y_p_zero y_pp_zero/2  y_ppp_zero/6  y_pppp_zero/24 0*y_5p_zero/factorial(5) 0*y_6p_zero/factorial(6)];

% coefficient = [c y_p_zero y_pp_zero/2  y_ppp_zero/6  y_pppp_zero/24];

% y_0_sol = simplify(coefficient(1) + y_pp_zero*s + y_pp_zero*s^2/2 + y_ppp_zero*s^3/6 + y_pppp_zero*s^4/24);

% f_1_sol = subs(y_p,y,y_0_sol);
               
% y_1_sol= c + int(f_1_sol,0,x);


if dece
   coefficient_num = double(subs(coefficient,[a b c gs],[miu_g slope max(Ux_backInt) g*sin(theta)]));

else
   coefficient_num = double(subs(coefficient,[a b c d gs ],[miu_g slope Ux_0 Ls g*sin(theta)]));

end


Ux_approximation = coefficient_num(1) + coefficient_num(2).*station_span + coefficient_num(3).*station_span.^2 + ...
                   coefficient_num(4) .*station_span.^3 + coefficient_num(5) .*station_span.^4 + ...
                   0*coefficient_num(6) .*station_span.^5 + 0*coefficient_num(7) .*station_span.^6;

               
               
h_fig = figure(212);
set(h_fig,'Name','path');
clf;
hold on
% plot(station_span,Ux_curceLimit,'b','LineWidth',2)
plot(station_span,Ux_frictionCircle_real,'g','LineWidth',2)

if dece
  plot(station_span,flip(Ux_approximation),'b--','LineWidth',2)

else
 plot(station_span,(Ux_approximation),'b--','LineWidth',2)

end

grid on
xlabel('station')
ylabel('velocity (m/s)')
% ylim([20 70])
% axis equal
legend('Numerical Solution','Polynomial Approximation')
box on
%% acceleration

% ax = (1./Ux_frictionCircle_real).* sqrt((miu_g.^2 - (slope*(Ls-station_span)'.*Ux_frictionCircle_real.^2).^2));
ay = slope*(Ls-station_span)'.*Ux_frictionCircle_real.^2;

ax = [sqrt(miu_g.^2 - ay(1).^2) ;0.5.*diff(Ux_frictionCircle_real.^2)./diff(station_span')];
% ay = v.^2.*path.curvature';

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
legend('ax','ay','a max')
% axis equal
box on


%% ========== section 3: given the paramers of a spiral, calculate the stopping distance  =======
%%

%% stop distance vs friction (acceleration from spiral end where curvature is 1/R)

% TRY backward integration
miu_g = (0.2:0.1:1)' * 9.81;
% miu_g = 0.8 * 9.8;
R_end = 200;

Ux_init = 0.1;
Ux_final = 35;
Ux = Ux_init:0.01:Ux_final;

% station_backInt = flip(station_spiral);
delta_Ux = diff(Ux);

S_stop = zeros(length(miu_g),length(Ux));

for i = 1:length(S_stop)-1
    S_stop(:,i+1) = S_stop(:,i) + (delta_Ux(i)*Ux(i))./ sqrt(miu_g.^2 - (slope.*(Ls-S_stop(:,i)).*Ux(i).^2).^2);
 
end

h_fig = figure(2100);
set(h_fig,'Name','path');
clf;
hold on
plot(Ux,S_stop,'b','LineWidth',2)
grid on
xlabel('Ux')
ylabel('Station (m)')
% ylim([20 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on

h_fig = figure(21001);
set(h_fig,'Name','path');
clf;
hold on
plot(miu_g/9.8,S_stop(:,end),'b','LineWidth',2)
grid on
xlabel('Friction Coef')
ylabel('Station (m)')
% ylim([20 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on

%% stop distance vs friction (deceleration from spiral start where curvature is 0)

% TRY backward integration
miu_g = (0.2:0.1:1)' * 9.81;
% miu_g = 0.8 * 9.8;
R_end = 200;

Ux_init = 35;
Ux_final = 0;
Ux = Ux_init:-0.01:Ux_final;

% station_backInt = flip(station_spiral);
delta_Ux = diff(Ux);

S_stop = zeros(length(miu_g),length(Ux));

for i = 1:length(S_stop)-1
    S_stop(:,i+1) = S_stop(:,i) - (delta_Ux(i)*Ux(i))./ sqrt(miu_g.^2 - (2*slope.*(S_stop(:,i)).*Ux(i).^2).^2);
 
end


h_fig = figure(2199);
set(h_fig,'Name','path');
clf;
hold on
plot(S_stop,Ux,'-','LineWidth',2)
grid on
ylabel('Ux')
xlabel('Station(m)')
% ylim([20 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on

h_fig = figure(21991);
set(h_fig,'Name','path');
clf;
hold on
plot(miu_g/9.8,S_stop(:,end),'b','LineWidth',2)
grid on
xlabel('Friction Coefficient')
ylabel('Stop Distance  (m)')
% ylim([20 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on

%% stop distance vs initial speed (deceleration from spiral start)

% TRY backward integration
% miu_g = (0.2:0.1:1)' * 9.81;
miu_g = 0.2 * 9.8;
R_end = 200;

clear Ux_init  Ux
Ux_init = 1:3:38;
Ux_final = 0;
for i = 1: length(Ux_init)
    Ux(i,:) = linspace(Ux_init(i),Ux_final,1000);
end
% Ux = Ux_init:-0.01:Ux_final;

% station_backInt = flip(station_spiral);
delta_Ux = diff(Ux,1,2);

S_stop = zeros(size(Ux));

for i = 1:length(S_stop)-1
    S_stop(:,i+1) = S_stop(:,i) - (delta_Ux(:,i).*Ux(:,i))./ sqrt(miu_g.^2 - (slope.*(S_stop(:,i)).*Ux(:,i).^2).^2);
 
end

h_fig = figure(3199);
set(h_fig,'Name','path');
clf;
hold on
plot(Ux_init,S_stop(:,end),'-','LineWidth',2)
grid on
xlabel('Ux0')
ylabel('stop distance(m)')
% ylim([20 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on

h_fig = figure(21991);
set(h_fig,'Name','path');
clf;
hold on
plot(S_stop,Ux,'b.','LineWidth',2)
grid on
ylabel('Ux')
xlabel('station(m)')
% ylim([20 60])
% axis equal
% legend('Curve limit speed','vehicle dynamics limit speed')
box on


%%
syms y(x)
f = y*sin(y);
% G = functionalDerivative(f,y)
diff(f,x)





%% test
syms y(x) a b x
y_p = y^4+ x^4;
% G = functionalDerivative(f,y)
y_pp = diff(y_p,x)


y_ppp = diff(y_pp,x)


y_pppp = diff(y_ppp,x)

%% Implicit differentiation
syms x y %Declaring symbilic variables
F(x,y) = x^2 + x*y + y^2 - 100 %Declaring implicit function
% Using Implicit Function Theorem
dy_dx = - diff(F,x)/diff(F,y)
