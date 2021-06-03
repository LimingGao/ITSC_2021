%%%%%%%%%%%%%%%%%%%%%  Function fcn_sectionRefQueryBySectionId %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%     Generate a Spiral curve given Ls, Rc, and delta_s.
%
% Input Variables:
%      Ls        Length measured along the spiral curve from its initial
%      position(curvature is zero) to end.(scalar)
%      Rc        Radius of circular curve at the end of the spiral.(scalar)
%      delta_s   The curve station resolution  
%
% Returned Results:
% 
% 
% Example:
% [x_spiral,y_spiral,station_spiral,slope] = fnc_spiralCurve(145,200);
%
% Processing Flow:
%
%
% Restrictions/Notes:
%
% The following functions are called:
%      none
%
% Author:             Liming Gao
% Created Date:       2021-04-01
% Revisions:
%     2021-06-02   add an input argument:delta_s
% 
% Reference: 
%     https://en.wikipedia.org/wiki/Euler_spiral
%
% To do list: nargin
% 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_spiral,y_spiral,station_spiral,slope] = fnc_spiralCurve(Ls,Rc,delta_s)

% DeNormallized Euler curve %https://en.wikipedia.org/wiki/Euler_spiral
% Ls = 45;%	Length measured along the spiral curve from its initial position
% Rc = 200; %Radius of circular curve at the end of the spiral

if nargin < 3 || isempty(delta_s) || ~isnumeric(delta_s)
    delta_s = 0.2; % the default resolution is 0.2 m
end

% Rc = 1/(2*Ls); %Radius of circular curve at the end of the spiral

fun_c = @(s) cos(s.^2);
fun_s = @(s) sin(s.^2);

factor_a= 1/sqrt(2*Rc*Ls); % a constant scaling up (denormalize) factor
slope = 2*factor_a^2;
% curvatre = 2*factor_a^2 * Ls
s_spiral = 0:delta_s*factor_a:Ls*factor_a; % scaled station

x_spiral = zeros(1, length(s_spiral));
y_spiral = zeros(1, length(s_spiral));
for n = 2:length(s_spiral)
    x_spiral(n) = x_spiral(n-1)+ integral(fun_c, s_spiral(n-1), s_spiral(n))/factor_a;
    y_spiral(n) = y_spiral(n-1)+ integral(fun_s, s_spiral(n-1), s_spiral(n))/factor_a;
    
end

station_spiral = [0 cumsum(sqrt(diff(x_spiral).^2+diff(y_spiral).^2))];

end


