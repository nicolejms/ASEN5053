% Patched Conics calculation script
clear
init
%% Parameters
% Constants
au = 1.496e8; % km
G = 6.6743e-20; % N*km^3/kg^2

% Masses
muSun = 1.3271e11; % km^3/s^2
muEarth = 3.986e5; % km^3/s^2
muJupiter = Gsaturn; % km^3/s^2
massEuropa = 480e20; % kg
muEuropa = G*massEuropa; % km^3/s^2

% Radii
REarth = 6371; % km
RJupiter = 69911; % km
REuropa = 1561; % km

% Orbital semi-major axes
aEarth = 1 * au; % km, around Sun
aJupiter = 3.2 * au; % km, around Sun
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

TOF = PeriodfromA(muSun,aHohEarthToJupiter)/2
Pdays = TOF/Day

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

% Hyperbolic eccentricities
eHypEarth = ecchfromArp(aHypEarth,rEarthSC)
eHypJupiter = EccfromArp(aHypJupiter,rJupiterSC)

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

% determine phase angles
    mmorigin = MMfromMua(muSun,aEarth);
    mmsc = MMfromMua(muSun,(aEarth+aJupiter)/2);      % mean motion of transfer 
    mmdestination = MMfromMua(muSun,aJupiter); % mean motion of destination
    
    phasei = pi -( mmdestination * TOF);  % phase angle for intercept 
    phasedegi = rad2deg(phasei);
    phasef = pi -( mmorigin * TOF);  % phase angle for intercept 
    phasedegf = rad2deg(phasef);


%% bielliptic transfer 

aIntermed = 2.3519e+08; 
% Earth-to-Intermediate orbit bielliptic transfer
aBiEarthToIntermed = (aEarth+aIntermed)/2;
vBiEarth = sqrt(2*muSun/aEarth - muSun/aBiEarthToIntermed);
vBiInterIn = sqrt(2*muSun/aIntermed - muSun/aBiEarthToIntermed);

% Intermediate-to-Jupiter bielliptic transfer
aBiIntermedToJupiter = (aJupiter+aIntermed)/2;
vBiIntermed = sqrt(2*muSun/aIntermed - muSun/aBiIntermedToJupiter);
vBiJupiterIn = sqrt(2*muSun/aJupiter - muSun/aBiIntermedToJupiter);

% Jupiter-to-Europa Hohmann transfer
aJupiterToEuropa = (rJupiterSC+aEuropa)/2;
vBiJupiterOut = sqrt(2*muJupiter/rJupiterSC - muJupiter/aJupiterToEuropa);
vBiEuropa = sqrt(2*muJupiter/aEuropa - muJupiter/aJupiterToEuropa);

% Hyperbolic excess velocities
vInfBiEarth = vBiEarth-vCircEarth;
vInfBiIntermed = vBiIntermed-vBiInterIn;
vInfBiJupiter = vCircJupiter-vBiJupiterIn;
vInfEuropa = vHohEuropa-vCircEuropa;

% Hyperbolic semi-major axes
aBiHypEarth = -muEarth/vInfBiEarth^2;
aBiHypJupiter = -muJupiter/vInfBiJupiter^2;
aHypEuropa = -muEuropa/vInfEuropa^2;

% Hyperbolic eccentricities
eHypBiEarth = ecchfromArp(aBiHypEarth,rEarthSC);
eHypBiJupiter = EccfromArp(aBiHypJupiter,rJupiterSC);

% Hyperbolic periapsis velocities
vHypBiPeriEarth = sqrt(2*muEarth/rEarthSC - muEarth/aBiHypEarth);
vHypBiPeriJupiter = sqrt(2*muJupiter/rJupiterSC - muJupiter/aBiHypJupiter);
vHypBiPeriEuropa = sqrt(2*muEuropa/rEuropaSC - muEuropa/aHypEuropa);

% Delta-v's
dvBiEarth = vHypBiPeriEarth-vCircEarthSC;
dvBiJupiter = abs(vBiJupiterOut-vHypBiPeriJupiter);
dvBiEuropa = abs(vCircEuropaSC-vHypBiPeriEuropa);

% Total Delta-v
dvBiTotal = dvBiEarth+vInfBiIntermed+dvBiJupiter+dvBiEuropa;
disp(dvBiTotal);

