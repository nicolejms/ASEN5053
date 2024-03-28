% Patched Conics calculation script

%% Parameters
% Constants
au = 1.496e8; % km
G = 6.6743e-20; % N*km^3/kg^2

% Masses
muSun = 1.3271e11; % km^3/s^2
muEarth = 3.986e5; % km^3/s^2
muJupiter = 1.2671e8; % km^3/s^2
massEuropa = 480e20; % kg
muEuropa = G*massEuropa; % km^3/s^2

% Radii
REarth = 6371; % km
RJupiter = 69911; % km
REuropa = 1561; % km

% Orbital semi-major axes
aEarth = 1 * au; % km, around Sun
aJupiter = 5.204 * au; % km, around Sun
aEuropa = 671000; % km, around Jupiter

% Spacecraft trajectory parameters
hEarth = 400; % km, LEO altitude
hEuropa = 400; % km, Europa orbit altitude
hJupiter = 1e6; % km, Jupiter flyby altitude
rEarthSC = REarth + hEarth; % km, spacecraft orbital radius
rJupiterSC = RJupiter + hJupiter; % km, spacecraft flyby radius
rEuropaSC = REuropa + hEuropa; % km, spacecraft orbital radius


%% Calculations
% Planetary circular orbit velocities
vCircEarth = sqrt(muSun/aEarth);
vCircJupiter = sqrt(muSun/aJupiter);
vCircEuropa = sqrt(muJupiter/aEuropa);

% Spacecraft circular orbit velocities
vCircEarthSC = sqrt(muEarth/rEarthSC);
vCircEuropaSC = sqrt(muEuropa/rEuropaSC);

% Earth-to-Jupiter Hohmann transfer
aHohEarthToJupiter = (aEarth+aJupiter)/2;
vHohEarth = sqrt(2*muSun/aEarth - muSun/aHohEarthToJupiter);
vHohJupiterIn = sqrt(2*muSun/aJupiter - muSun/aHohEarthToJupiter);

% Jupiter-to-Europa Hohmann transfer
aJupiterToEuropa = (rJupiterSC+aEuropa)/2;
vHohJupiterOut = sqrt(2*muJupiter/rJupiterSC - muJupiter/aJupiterToEuropa);
vHohEuropa = sqrt(2*muJupiter/aEuropa - muJupiter/aJupiterToEuropa);

% Hyperbolic excess velocities
vInfEarth = vHohEarth-vCircEarth;
vInfJupiter = vCircJupiter-vHohJupiterIn;
vInfEuropa = vHohEuropa-vCircEuropa;

% Hyperbolic semi-major axes
aHypEarth = -muEarth/vInfEarth^2;
aHypJupiter = -muJupiter/vInfJupiter^2;
aHypEuropa = -muEuropa/vInfEuropa^2;

% Hyperbolic periapsis velocities
vHypPeriEarth = sqrt(2*muEarth/rEarthSC - muEarth/aHypEarth);
vHypPeriJupiter = sqrt(2*muJupiter/rJupiterSC - muJupiter/aHypJupiter);
vHypPeriEuropa = sqrt(2*muEuropa/rEuropaSC - muEuropa/aHypEuropa);

% Delta-v's
dvEarth = vHypPeriEarth-vCircEarthSC;
dvJupiter = abs(vHohJupiterOut-vHypPeriJupiter);
dvEuropa = abs(vCircEuropaSC-vHypPeriEuropa);

% Total Delta-v
dvTotal = dvEarth+dvJupiter+dvEuropa;
disp(dvTotal);




