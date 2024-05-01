function rrd = jpleph_mice (et, ntarg, ncent)

% reads the jpl planetary ephemeris and gives the position and velocity
% of the point 'ntarg' with respect to point 'ncent' using MICE routines

% input

%   et    = TDB julian date at which interpolation is wanted

%   ntarg = integer number of 'target' point

%   ncent = integer number of center point

%   the numbering convention for 'ntarg' and 'ncent' is:

%        1 = mercury (199)     8 = neptune BC (8)
%        2 = venus (299)       9 = pluto BC (9)
%        3 = earth (399)      10 = moon (301)
%        4 = mars BC (4)      11 = sun (10)
%        5 = jupiter BC (5)   12 = solar system barycenter (0)
%        6 = saturn BC (6)
%        7 = uranus BC (7)    BC ==> barycenter

% output

%   rrd = output 6-word array containing position and velocity
%         of point 'ntarg' relative to 'ncent'. the units are
%         determined by the value of km passed via global.

% global

%   iephem  = initialization flag (1 = initialize)
%   ephname = name of ephemeris binary data file (de421.bsp, etc.)
%   km      = state vector units flag 
%             (1 = kilometers and kilometers/second)
%             (0 = aunits and aunits/day)
%   aunit   = numerical value of astronomical unit (kilometers)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iephem ephname km aunit

if (iephem == 1)
    
    % load binary ephemeris data file
    
    cspice_furnsh(ephname);
    
    % reset initialization flag
    
    iephem = 0;
    
end

% set name of target body

switch ntarg
    
    case (1)
        
        targ = 'mercury';
        
    case (2)
        
        targ = 'venus';
        
    case (3)
        
        targ = 'earth';
        
    case (4)
        
        % Mars barycenter
        
        targ = '4';
        
    case(5)
        
        % Jupiter barycenter
        
        targ = '5';
        
    case (6)
        
        % Saturn barycenter
        
        targ = '6';
        
    case (7)
        
        % Uranus barycenter
        
        targ = '7';
        
    case (8)
        
        % Neptune barycenter
        
        targ = '8';
        
    case (9)
        
        % Pluto barycenter
        
        targ = '9';
        
    case (10)
        
        targ = 'moon';
        
    case (11)
        
        targ = 'sun';
        
    case (12)
        
        % solar system barycenter
        
        targ = '0';
        
end

% set name of central body

switch ncent
    
    case (1)
        
        obs = 'mercury';
        
    case (2)
        
        obs = 'venus';
        
    case (3)
        
        obs = 'earth';
        
    case (4)
        
        % Mars barycenter
        
        obs = '4';
        
    case (5)
        
        % Jupiter barycenter
        
        obs = '5';
        
    case (6)
        
        % Saturn barycenter
        
        obs = '6';
        
    case (7)
        
        % Uranus barycenter
        
        obs = '7';
        
    case (8)
        
        % Neptune barycenter
        
        obs = '8';
        
    case (9)
        
        % Pluto barycenter
        
        obs = '9';
        
    case (10)
        
        obs = 'moon';
        
    case (11)
        
        obs = 'sun';
        
    case (12)
        
        % solar system barycenter
        
        obs = '0';
        
end

% compute time, expressed as TDB seconds past J2000 TDB (2451545.0)

etime = 86400.0d0 * (et - 2451545.0d0);

% compute position and velocity vectors in eme2000 ecliptic system (no corrections)

starg = mice_spkezr(targ, etime, 'ECLIPJ2000', 'NONE', obs);

% provide output in user-requested units

if (km == 1)
    
    % state is kilometers and kilometers/second
    
    rrd = starg.state;
    
else
    
    % state is au's and au's/day
    
    rrd(1:3) = starg.state(1:3) / aunit;
    
    rrd(4:6) = 86400.0 * starg.state(4:6) / aunit;
    
    rrd = rrd';
    
end





