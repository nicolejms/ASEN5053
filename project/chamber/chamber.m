clc;
clear all;
close all;

g0 = 9.81; %m/s
m_dot = 78.6659; %mass flow rate of propellant [kg/s]
max_stress = 2.48211e8; %allowable working stress for steel in Pascals
safety_factor = 2;
rho_wall = 7850; %Density of the wall material assuming steel [kg/m^3]

Contraction_ratio = 3; %Input desired contraction ratio Ac/At (usually between 3-5)
Contraction_angle = 45; %Degrees
cone_half_angle = 15; %Degrees

L_star = 0.762; %Characteristic length in meters. NTO/MMH has a range of 0.762-0.889[m]

Index_P = 8; %Desired index for pressure (8 for 80 atm)
Index_expansion_ratio = 3;

Pc = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]; %Chamber pressure range [atm]
Pc_pascal = Pc*101325; %Chamber pressure in Pascals
c_vac = [2140.3, 2153.9, 2161.3, 2166.3, 2170.0, 2173.0, 2175.4, 2177.4, 2179.4, 2180.6, 2182.0, 2183.2, 2184.3, 2185.2, 2186.2];
Isp_vac = c_vac./g0;

%Expansion Ratio Range for O/F = 1.9 and Pc = 80 atm
expansion_ratio_range = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130];
c_vac_expanded = [2973.2,3090.3,3146.0,3180.8,3205.3,3223.9,3238.7,3250.8,3261.1,3269.9,3277.5,3284.2,3290.3];
Isp_vac_expanded = c_vac_expanded/g0; %Isp_vac in  seconds


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
m_CC = pi.*rho_wall.*actual_thickness_CC.*(2.*r_cc(Index_P).*L_CC + (pi.*(r_cc(Index_P).^2 - r_t(Index_P).^2)/(tand(Contraction_angle))));

%Calculating Nozzle Length and Mass (assuming wall thickness is the same as
%combustion chamber and a conical nozzle)
A_exit = expansion_ratio_range.*A_throat(Index_P);
diameter_exit = sqrt((4.*A_exit)/pi); %Nozzle exit diameter [m]
r_exit = diameter_exit./2; %Nozzle exit radius [m]
L_nozzle = (r_exit - r_t(Index_P))/tand(cone_half_angle); %Nozzle length [m]

m_nozzle = pi*rho_wall*actual_thickness_CC(Index_P).*L_nozzle.*(r_exit + r_t(Index_P));

%Chamber Dimensions based on 80 atm chamber pressure
fprintf('Chamber Pressure = %.2f atm\n', Pc(Index_P));
fprintf('Throat Area = %.5f m^2\n', A_throat(Index_P));
fprintf('Throat Diameter = %.5f m\n', diameter_throat(Index_P));
fprintf('Combustion Chamber Area = %.5f m^2\n', A_CC(Index_P));
fprintf('Combustion Chamber Diameter = %.5f m\n', diameter_CC(Index_P));
fprintf('Combustion Chamber Volume = %.5f m^3\n', V_cc(Index_P));
fprintf('Combustion Chamber Length = %.5f m\n', L_CC);
fprintf('Chamber Wall Thickness = %.5f m\n', actual_thickness_CC(Index_P));
fprintf('Chamber Mass = %.5f kg\n', m_CC(Index_P));
%Nozzle Dimensions based on expansion ratio
fprintf('Expansion Ratio = %.2f\n', expansion_ratio_range(Index_expansion_ratio));
fprintf('Nozzle Exit Area = %.5f\n', A_exit(Index_expansion_ratio));
fprintf('Nozzle Exit Diameter = %.5f m\n', diameter_exit(Index_expansion_ratio));
fprintf('Nozzle Length = %.5f m\n', L_nozzle(Index_expansion_ratio));
fprintf('Nozzle Mass = %.5f kg\n', m_nozzle(Index_expansion_ratio));
fprintf('Isp_{vac} = %.5f sec\n', Isp_vac_expanded(Index_expansion_ratio));
fprintf('Effective Exhaust Velocity (Vacuum) = %.2f m/s\n', c_vac_expanded(Index_expansion_ratio));

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

% Plotting Isp vs Expansion Ratio
figure;
yyaxis left; % Set the left y-axis for specific impulse
plot(expansion_ratio_range, Isp_vac_expanded);
xlabel('Expansion Ratio');
ylabel('Isp_{vac} [sec]');
title('Isp_{vac} of NTO/MMH vs. Expansion Ratio');
grid on;

hold on;

yyaxis right; % Set the right y-axis for wall thickness
plot(expansion_ratio_range, m_nozzle, '--');
ylabel('Nozzle Wall Mass [kg]');

legend('Isp_{vac}', 'Wall Mass','Location','southeast'); % Adding legend