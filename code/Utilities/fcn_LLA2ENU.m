%%%%%%%%%%%%%%%%%%%%%%%%%%% Function fcn_LLA2ENU %%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   fcn_LLA2ENU converts an array of points from LLA to ENU
% 
% Format:
%   enu = fcn_LLA2ENU(arrayLat, arrayLon, arrayAlt, refLat, refLon, refAlt)
% 
% INPUTS:
%   arrayLat, arrayLon, arrayAlt: array of points in LLA
%   refLat, refLon, refAlt: origin of the ENU coordinate system
% 
% OUTPUTS:
%   enu: 'arrayLat, arrayLon, arrayAlt' in ENU coordinate system
% 
% NOTE: Requires function rot to be in the same directory
% Author: Satya Prasad
% Created: 2020-06-20
% 
% Reference:
% https://gssc.esa.int/navipedia/index.php/Reference_Frames_in_GNSS
% https://gssc.esa.int/navipedia/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion
% https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
% 
% Revision history:
% 
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function enu = fcn_LLA2ENU(arrayLat, arrayLon, arrayAlt, refLat, refLon, refAlt)
%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are there right number of inputs?
if nargin ~= 6
    error('fcn_LLA2ENU: Incorrect number of input arguments');
end

if ((refLat < -90.0) || (refLat > +90.0) || (refLon < -180.0) || (refLon > +180.0))
    error('fcn_LLA2ENU: WGS reflat or WGS reflon - out of range');
end

if numel(refLat) ~= 1
    error('fcn_LLA2ENU: reflat must be a scalar');
end

if numel(refLon) ~= 1
    error('fcn_LLA2ENU: reflon must be a scalar');
end

if numel(refAlt) ~= 1
    error('fcn_LLA2ENU: refalt must be a scalar');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants
D0         = pi/180;            % conversion deg to rad
A_EARTH    = 6378137;           % Semi-major axis of the ellipse in meters
flattening = 1/298.257223563;   % flattening factor
% square of the first numerical eccentricity of the ellipsoid
NAV_E2 = (2-flattening)*flattening;

lat = [refLat; arrayLat];
lon = [refLon; arrayLon];
alt = [refAlt; arrayAlt];
% compute sine and cosine of latitude
slat = sin(D0*lat);
clat = cos(D0*lat);
% compute sine and cosine of longitude
slon = sin(D0*lon);
clon = cos(D0*lon);

% radius of curvature in the prime vertical
r_n = A_EARTH./sqrt(1 - NAV_E2*slat.^2);

% ECEF coordinate system
xyz = [ (r_n + alt) .* clat .* clon, ...
        (r_n + alt) .* clat .* slon, ...
        ((1 - NAV_E2)*r_n + alt) .* slat ];

% reflat, reflon, refalt in ECEF coordinate system
refxyz = xyz(1,:);

% Vector pointing from reference point to the target in ECEF coord system
diffxyz = xyz(2:end,:) - refxyz;

% Now rotate the diffxyz vector to enu frame
R1 = rot(90-refLat, 1); % Rotation about e1 axis
R3 = rot(90+refLon, 3); % Rotation about e3 axis
R = R1*R3;

enu = [ sum(R(1,:).*diffxyz, 2), ...
        sum(R(2,:).*diffxyz, 2), ...
        sum(R(3,:).*diffxyz, 2) ];

end