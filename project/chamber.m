clc;
%clear all;
close all;

g0 = 9.81; %m/s
m_dot = 78.6659; %mass flow rate of propellant [kg/s]
max_stress = 2.48211e8; %allowable working stress for steel in Pascals
safety_factor = 2;
rho_wall = 7850; %Density of the wall material assuming steel [kg/m^3]

Contraction_ratio = 3; %Input desired contraction ratio Ac/At (usually between 3-5)
Contraction_angle = 45; %Degrees

L_star = 0.762; %Characteristic length in meters. NTO/MMH has a range of 0.762-0.889[m]

Pc = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]; %Chamber pressure range [atm]
Pc_pascal = Pc*101325; %Chamber pressure in Pascals
c_vac = [2140.3, 2153.9, 2161.3, 2166.3, 2170.0, 2173.0, 2175.4, 2177.4, 2179.4, 2180.6, 2182.0, 2183.2, 2184.3, 2185.2, 2186.2];
Isp_vac = c_vac./g0;

%Chamber Parameters based on an O/F of 1.9
gamma = 1.1558; %Specific heat ratio in the chamber for an O/F of 1.9
Tc = 3344.26; %Chamber temperature
M = .02186; %Molecular mass in kg/mol
R = 8.3145/M; %Specific gas constant in J/kg-K

%Throat area and diameter based on mass flow rate of propellant
A_throat = (m_dot./(Pc_pascal.*gamma))*((sqrt(gamma*R*Tc))/(sqrt(((2/(gamma+1))^((gamma+1)/(gamma-1)))))); %[m^2]
diameter_throat = sqrt((4.*A_throat)/pi); %Throat diameter [m]
r_t = diameter_throat./2; %Throat radius [m]

%Combustion Chamber Length and Area
A_CC = A_throat.*Contraction_ratio; %[m^2]
diameter_CC = sqrt((4.*A_CC)/pi); %Combustion Chamber diameter [m]
r_cc = diameter_CC./2; %Combustion Chamber radius [m]
L_CC = L_star.*(A_throat/A_CC); %Combustion Chamber length in [m]

%Volume of combustion chamber
V_cc = L_star.*A_throat; %[m^3]

%Thickness
thickness_CC = (Pc_pascal.*diameter_CC)/(2*max_stress);
actual_thickness_CC = thickness_CC*safety_factor; %Thickness in [m]

%Mass of Combustion Chamber based on largest volume calculation
m_CC = pi.*rho_wall.*actual_thickness_CC.*(2.*r_cc(1).*L_CC + (pi.*(r_cc(1).^2 - r_t(1).^2)/(tan(Contraction_angle))));

% Plotting Isp vs Chamber Pressure
figure;
yyaxis left; % Set the left y-axis for specific impulse
plot(Pc, Isp_vac);
xlabel('Chamber Pressure [atm]');
ylabel('Isp_{vac} [sec]');
title('Isp of NTO/MMH vs. Chamber Pressure');
grid on;

hold on;

yyaxis right; % Set the right y-axis for wall thickness
plot(Pc, actual_thickness_CC, '--');
ylabel('Wall Thickness [m]');

legend('Isp_{vac}', 'Wall Thickness','Location','southeast'); % Adding legend


% Plotting Chamber Volume vs Chamber Pressure
figure;
yyaxis left; % Set the left y-axis for specific impulse
plot(Pc, Isp_vac);
xlabel('Chamber Pressure [atm]');
ylabel('Isp_{vac} [sec]');
title('Isp of NTO/MMH vs. Chamber Pressure');
grid on;

hold on;

yyaxis right; % Set the right y-axis for wall thickness
plot(Pc, m_CC, '--');
ylabel('Combustion Chamber Wall Mass [kg]');

legend('Isp_{vac}', 'Wall Mass','Location','southeast'); % Adding legend