clc
%clear
close all 

% mission parameters
OFratio = 1.9;          % oxidizer to fuel ratio
Ftotal = 252577;        % total thrust N
RatioContraction = 3;   % contraction ratio between chamber and throat area
angle_conv = 45;        % angle of convergence
Isp = 325;              % specific impulse sec
TstayMin = 0.01;        % min 10 ms time in chamber

% constants
Rc = 8.3144;            % gas contant J/ Kmol
g0 = 9.806;             % gravitational constant m/s^2
MwMmh = 0.046074 ;      % mol weight mmh
MwNto = 0.09201;        % mol weight nto
k_mmh = 1.1558;         % ratio of specific heats for MM
atm2pascal = 101325;    % atm to pascal conversion factor

% from CEARun
k_mixture = 1.1537;     % from CEARun
MwMix = 0.02186;        % kg/mol
R = Rc/MwMix            % J/kg K 

P_chamber_atm = 20:10:100;    % chamber pressure range 20-100 atm
P_chamber = P_chamber_atm.*atm2pascal;
T_chamber = [3220.22 3258.05 3284.18 3303.99 3319.86 3333.03 3344.26 3354.01 3362.61];
rho_chamber = [1.6366 2.4345 3.2275 4.0172 4.8043 5.5893 6.3726 7.1545 7.9351 ];
C_vac = []/g0;
% equations %

% Area of the throat
AtfromK = @(mdot,k,T1,P1,R) (mdot./(P1*k)).*((sqrt(k*R.*T1))/(sqrt(((2/(k+1))^((k+1)/(k-1)))))); %[m^2]
AtfromCstar = @(mdot,Cstar,P1) mdot*Cstar/P1;

% characteristic chamber length
LstarfromVol = @(Vol,At) Vol./At;
LstarfromMdot = @(mdot,At,rho,tstay) (mdot*tstay/g0)./(At.*rho);

% stay time
TstayfromMdot = @(Vol,mdot,density) Vol*density/mdot;

% chamber volume
VolCcfromAreas = @(A1,At,L1,Lc) A1*L1 + A1*Lc*(1+sqrt(At/A1) + At/A1);
VolCcfromTs = @(mdot,rho,tstay) mdot*tstay./(rho*g0);

% Thrust
FfromCstar = @(mdot,Cstar,Cf) mdot*Cstar*Cf;
FfromV2 = @(mdot,V2,P2,A2) mdot*V2 + P2*A2;

% Effective exhaust velocity
CfromV2 = @(V2,P2,mdot) V2 + P2/mdot;
CfromIsp = @(Isp) Isp *g0;

% mass flow rate
MdotfromAt = @(Pc,At,Cstar) Pc*At/Cstar;
MdotfromCF = @(F,C) F/C;

% characteristic velocity
CstarfromC = @(c,Cf) c/Cf;
CstarfromAt = @(P1,At,mdot) P1.*At/mdot;
CstarfromK = @(k,T1,R) sqrt(k*R*T1)/(k*sqrt((2/(k+1))^(k+1/(k-1))));

% calculations %

% find effective exhaust velocity
C = CfromIsp(Isp);

% find mass flow rate
mdot = Ftotal/C;

% calc area of throat and combustion chamber from k and mdot
At = AtfromK(mdot,k_mixture,T_chamber,P_chamber,R);  % area of throat
Dt = sqrt((4*At)/pi)    ;                            % diameter of throat
Acc = At * RatioContraction;                         % area of combustion chamber
Dcc = sqrt((4*Acc)/pi);                              % diameter of combustion chamber

% calc Cstar and Coeff of Thrust Cf
Cstar = CstarfromAt(P_chamber,At,mdot);
Cf = C./Cstar;

% calculate volume of combustion chamber and L1,Lstar based on stay time
VolCc = VolCcfromTs(mdot,rho_chamber,TstayMin);
Lstar = LstarfromVol(VolCc,At);
Lchamber = (VolCc./At - (sqrt(At/pi)*cotd(Contraction_angle)*(Contraction_ratio^(1/3)-1)/3))/Contraction_ratio;




% Plotting At vs Chamber Pressure
figure;
yyaxis left; % Set the left y-axis for specific impulse
plot(P_chamber_atm, At);
xlabel('Chamber Pressure [atm]');
ylabel('A_{throat} [m]');
title('Throat Area and Characteric Length for NTO/MMH vs. Chamber Pressure');
grid on;

hold on;

yyaxis right; % Set the right y-axis for wall thickness
plot(P_chamber_atm, Lstar, '--');
ylabel('L^* [m]');

legend('A_{throat}', 'L^*','Location','east'); % Adding legend

% Plotting Chamber Volume vs Chamber Pressure
figure;
yyaxis left; % Set the left y-axis for specific impulse
plot(P_chamber_atm, Cstar);
xlabel('Chamber Pressure [atm]');
ylabel('C^* [m]');
title('C^* and Thrust Coefficient of Chamber for NTO/MMH vs. Chamber Pressure');
grid on;

hold on;

yyaxis right; % Set the right y-axis for wall thickness
plot(P_chamber_atm, Cf, '--');
ylabel('Coefficient of Thrust - C_F');

legend('C^*', 'C_F','Location','east'); % Adding legend