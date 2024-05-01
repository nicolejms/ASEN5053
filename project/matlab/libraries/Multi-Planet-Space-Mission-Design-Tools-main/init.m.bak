%% rotation matrices
M1 = @(a) [1 0 0; 0 cosd(a) sind(a); 0 -sind(a) cosd(a)];
M2 = @(a) [];
M3 = @(a) [cosd(a) sind(a) 0; -sind(a) cosd(a) 0; 0 0 1];
% DCM  rth to xyz
DMC = @(Raan, i, theta) M3(-Raan) * M1(-i) * M3(-theta);

%% constants
x_hat = [1;0;0];
y_hat = [0;1;0];
z_hat = [0;0;1];

Minute = 60;
Hour = 60*60;
Day = 24*60*60;
Year = 365*24*60*60;

gc = 32.174;  % conversion
g0 = 32.174;  % accel due to gravity
psi2lbf = 144;  

% mass
Meuropa = 4.79984e22; 

% gravitational parameter
G = 6.673e-20;
Gmars = 4.305e4;
Gmoon = 4902.799;
Gsun = 1.32712428e11;
Gsaturn = 3.794e7;
Gjupiter = 1.268e8;
Gearth = 3.986004415e5;
Gvenus = 3.257e5;
Geuropa = G * Meuropa;

% radius
Rearth = 6378;  
Rmars = 3397.2;
Rmoon = 1738;
Rsaturn = 60268;
Rjupiter = 71492;
Rvenus = 6052.0;
Reuropa = 1565;

% distance from sun
AofEarthHelio = 1.0000010178;
AofSaturnHelio = 9.554909595;
AofJupiterHelio = 5.202603191;
AofMarsHelio = 1.524;   
AofVenus = 1.08208601e8; % km
AU = 149597870.7;

% orbit of Europa
periEuropa = 664862;    % km
apoEuropa = 676938;     % km
radEuropa = 670900;       % km

% sphere of influence 
SOIjupiter = 4.82e7;
SOIEarth = 9.25e5;
SOIMars = 5.78e5;
SaturnJ2 = 0.01645;

%% radius 
% radius of periapsis
Rp = @(a,e) a * (1 - e);

% radius of apoapsis
Ra = @(a,e) a * (1 + e);

% semi latus rectum
PfromAe = @(a,e) a*(1-e^2);

% radius
Rradial = @(mu,e,h,ta) (h^2/mu)*(1+e*cosd(ta)); % from ang momentum, true anomaly, ecc
RfromP = @(p,e,ta) p/(1+e*cosd(ta));        % from semi-latus rectum
RfromEA = @(a,e,ea) a*(1-e*cos(ea));        % from eccentric anomaly

% r vector - |r| (cos(raan) + sin(raan))
RfromAngle = @(r, angle) r*[cosd(angle); sind(angle); 0];
RfromConic = @(mu,h,e,ta) h^2/(1+e*cosd(ta))/mu;
RfromVEs = @(mu,v,E) mu/((v^2/2)-E);        % from velocity and specific energyy

%% velocity
VfromEs = @(mu,Es,r) sqrt(2*(Es+mu/r));     % from specific energy
Vradial = @(mu,e,h,ta) (mu/h) * e*sind(ta); % radial velocity
Vtrans = @(mu,e,h,ta) (mu/h)*(1+e*cosd(ta)); % transverse velocity
VfromRA = @(mu,r,a) sqrt((2*mu/r)-(mu/a));  % from radius and sma
Vescape = @(mu,r) sqrt(2*mu/r);             % escape velocity

%% orbital elements

% semi major axis
AfromVR = @(mu,v,r) -mu/(v^2 -(2*mu/r));    % from velocity, radius
AfromE = @(mu,E) -mu/(2*E);                 % from specific energy
AfromPeriod = @(mu,period) nthroot((period/(2*pi))^2 * mu,3);
AfromRp = @(periapsis,e) periapsis/(1-e);   % from radius of periapsis

% semi minor axis
BfromAE = @(a,e) a * sqrt(1-e^2);   

% eccentricity
EccfromHE = @(mu,h,E) sqrt(1+(2*h^2*E/mu^2));       % from specific energy and ang mom
EccfromHRV = @(mu,r,v,H) cross(v,H)/mu - r/norm(r); % from ang mom, velocity
EccfromArp = @(a,rp) 1-rp/a;                        % from rad of periapsis
EccfromAra = @(a,ra) ra/a-1;                        % from rad of apoapsis
EccfromRpRa = @(rp,ra) (ra-rp)/(ra+rp);
EccfromAP = @(a,p) sqrt(1-p/a);                     % from semi-latus rectum

% specific angular momentum
HfromRV = @(r,v) cross(r,v);
HfromEe = @(mu,e,E) sqrt((e^2-1)*mu^2/(2*E));
HfromP = @(mu,p) sqrt(mu*p);
HfromTa = @(mu,r,e,ta)sqrt(r*(1+e*cosd(ta))*mu);
HfromRe = @(mu,r,e) sqrt(mu*(norm(r) + dot(r, e)));
HfromAe = @(mu,a,e) sqrt(mu*a*(1-e^2));

% specific energy
EfromA = @(mu,a) -mu/(2*a); % energy from sma
EfromHe = @(mu,e,h) (e^2-1)*mu^2 / (2* h^2); % energy from angular momentum
EfromVr = @(mu,v,r) v^2/2 - (mu/r); % energy from velocity, position

% period
PeriodfromA = @(mu,a) 2*pi*sqrt(a^3/mu); % 

% inclination
IfromH = @(h) acosd(h(3)/norm(h));

% RAAN
LONfromH = @(h) cross(z_hat, h);
RaanfromLON = @(lon) acosd(lon(1)/norm(lon));
raan_sign = @(lon) sign(dot(lon, y_hat));

% argument of periapasis
WfromLONe = @(lon,e) acosd(dot(lon,e)/(norm(lon)*norm(e)));
% w sign check
w_sign = @(e) sign(dot(e,z_hat));

% true anomoly
TAfromRe = @(e,r) acosd(dot(e,r)/(norm(r)*norm(e)));
TAfromHr = @(mu,H,r,e) acosd((H^2/(mu*r) -1)/e);
TAfromAe = @(atransfer,aplanet,e) acosd(( ((atransfer/aplanet)*(1-e^2))-1)/e);
% true anomaly from Ecc Anomaly
TAfromEA = @(e,ea) atand((sqrt((1+e)/(1-e))*tan(ea/2)))*2;

% true anomaly from velocity
% TAfromVrad = @(mu,e,h,ta) (mu/h) * e*sind(ta);
% TAfromVtrans = @(mu,e,h,ta) (mu/h)*(1+e*cosd(ta)); 

% true anomoly sign check
ta_sign = @(r,v) sign(dot(r,v));

% eccentric anomaly
EAfromRae = @(r,a,e) acos((r/a-1)/(-e));
EAfromEta = @(e,ta) acos((e+cosd(ta))/(1+e*cosd(ta)));
EAfromTatan = @(e,ta) atan((sqrt((1-e)/(1+e))*tand(ta/2)))*2;

% Mean motion, n
MMfromMua = @(mu,a) sqrt(mu/(a^3));
MMfromP = @(p) 2*pi/p;

%% time in orbit
T_afromb = @(n,e,Ea,Eb) ((Eb-e*sin(Eb)) - (Ea-e*sin(Ea)))/n;
T_bfroma = @(n,e,Ea,Eb) ((Ea-e*sin(Ea)) - (Eb-e*sin(Eb)))/n;
T_pastPeri = @(n,e,E) (1/n)*(E-e*sin(E));

%% Impulse
Mfinal = @(Minit,deltaV,Is) Minit * exp(-(deltaV*1000)/(Is*9.81));
Minitial = @(Mfin,deltaV,Is) Mfin * exp((deltaV*1000)/(Is*9.81));

%% transfer orbits
% Hohmann

v1_init = @(mu,r1t) sqrt(mu/r1t);
v2_final = @(mu,r2f) sqrt(mu/r2f);
synodFromMm = @(mm1,mm2) abs(2*pi/(mm1-mm2));
synodFromMua = @(mu,a1,a2) 2*pi / (sqrt(mu/a1^3) - sqrt(mu/a2^3));

% interplanetary
TAfromPRE = @(p,r,e) acosd((1/e*((p/r)-1)));
vout = @(vinf,vplanet,angle) sqrt(vinf^2 + vplanet^2 - (2*vinf*vplanet*cosd(angle)));
deltavequiv = @(vinf, turningangle) 2*vinf*sind(turningangle/2);

EccfromParkOrbitVinf = @(mu,rpark,vinf) 1 + rpark * vinf^2 / mu;
burnAngle = @(ecc) acos(-1/ecc);

% parabolic transfer orbit
TOFp_plus = @(mu,s,c) 1/3 * sqrt(2/Gsun) * (s^(3/2)-(s-c)^(3/2));
TOFp_minus = @(mu,s,c) 1/3 * sqrt(2/Gsun) * (s^(3/2)+(s-c)^(3/2));

% f & g
F = @(r,p,theta) 1-(r/p) * (1-cosd(theta));
G = @(mu,r,r0,p,theta) r*r0 * sind(theta) /sqrt(mu*p);
Fdot = @(mu,r,r0,p,theta) sqrt(mu/p)*tand(theta/2)*((1-cosd(theta))/p - 1/r - 1/r0);
Gdot = @(r,p,theta) 1-(r/p)*(1-cosd(theta));

fFromdeltaEA = @(sma,r0,EA) 1 - sma/r0 * (1-cos(EA));
gFromdeltaEA = @(mu,sma,r0,EA,dt) dt - sqrt((sma^3/mu)) * (EA-sin(EA));
fdotFromdeltaEA = @(mu,sma,r0,r,EA) (-sin(EA) * sqrt(mu*sma))/(r0*r);
gdotFromdeltaEA = @(sma,r,EA) 1 - (sma/r) * (1-cos(EA));

% hyperbolas
ecchfromArp = @(a,rp) 1-(rp/a);
delta = @(e) 2*asind(1/e);

% perturbations
muFromPa = @(period,a)  (4*pi^2*a^3)/(period^2);
J2 = @(mu,a,e,i,R,delRaan) (-2*delRaan*(1-e^2)^2*a^(7/2))/(3*sqrt(mu)*R^2*cosd(i));
iFromJ2 = @(mu,a,e,R,delRaan,J2) acosd((delRaan*(-2/3)*(1-e^2)^2*a^(7/2))/(sqrt(mu)*J2*R^2));


% relative trajectories and Hill equations
radial_x = @(x0,y0,z0,xdot0,ydot0,zdot0,n,t) 4*x0 +2*ydot0/n + xdot0*sin(n*t)/n - (3*x0+2*ydot0/n)*cos(n*t);
along_y = @(x0,y0,z0,xdot0,ydot0,zdot0,n,t) y0 -2*xdot0/n -3*(3*x0 +2*ydot0/n)*sin(n*t) + 2*xdot0*cos(n*t)/n;
 cross_z = @(z0,zdot0,n,t) zdot0*sin(n*t)/n + z0*cos(n*t);

rendezvous_rr = @(n,t) [4-3*cos(n*t) 0 0 ;
    6*(sin(n*t)-n*t) 1 0;
    0 0 cos(n*t)];

rendezvous_rv = @(n,t) [sin(n*t)/n 2*(1-cos(n*t))/n 0;
    2*(cos(n*t)-1)/n 4*sin(n*t)/n-3*t 0;
    0 0 sin(n*t)/n];

rendezvous_vr = @(n,t) [3*n*sin(n*t) 0 0;
    6*n*(cos(n*t)-1) 0 0;
    0 0 -n*sin(n*t)];

rendezvous_vv = @(n,t) [cos(n*t) 2*sin(n*t) 0;
    -2*sin(n*t) 4*cos(n*t)-3 0;
    0 0 cos(n*t)];
