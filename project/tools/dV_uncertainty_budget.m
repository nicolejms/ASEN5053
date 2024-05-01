%% Propellant Budget

% Requirement based on required dV

% Key Criteria (Assume worst case scenario)
% Main Flight Chamber: 85% - 96%
% Flight Control: 2% - 10%
% Residual Propellant: 0.5% - 2%
% Loading Uncertainty: up to 0.5%
% Allowance for off-nominal performance: 0.1% - 1%
% Allowance for off-nominal operations: 0.1% - 1%
% Mission Margin: 1% - 5% of Main Flight Chamber & Flight Control
% Contingency: 1% - 5% of total load

dV = 14.1979; % km/s (Assumed 85% of total DV budget, main flight chamber)
dV_flight_tot = (1/0.85)*dV; % km/s (total DV needed for all flight operations)
dV_FC = 0.10*dV_flight_tot; % km/s (dV needed for flight control, assumed 10%)
dV_RP = 0.02*dV_flight_tot; % km/s (dV needed to account for residual, unused propellant, assumed 2%)
dV_LU = 0.005*dV_flight_tot; % km/s (dV to account for loading uncertainty, assumed 0.5%)
dV_ONP = 0.01*dV_flight_tot; % km/s (dV to account for off nominal performance, assumed 1%)
dV_ONO = 0.01*dV_flight_tot; % km/s (dV to account for off nominal operations, assumed 1%)

dV_MM = (dV+dV_FC)*0.05; % km/s (dV alloted for mission margins, assumed 5%)
dV_cont = dV_flight_tot*0.05; % km/s (dV alloted for total mission contingency, assumed 5%)
dV_max = dV_flight_tot + dV_cont; % km/s (Maximum dV requirement)

FM = 1000; % Final spacecraft mass, kg
c = 3.1889; % Effective exhaust velocity, km/s
dV_vec = [dV, dV_flight_tot, dV_FC, dV_RP, dV_LU, dV_ONP, dV_ONO, dV_MM, dV_cont, dV_max]; % Vector of allocated dVs
Prop_mass_allo = FM.*(exp(dV_vec./c)-1); % Vector of propellant mass allocations, kg

OFR = 2; % Oxidizer to fuel ratio
Prop_mass_ox = (OFR/(OFR+1)).*Prop_mass_allo; % Vector of oxidizer mass allocations, kg
Prop_mass_fuel = (1/(OFR+1)).*Prop_mass_allo; % Vector of fuel mass allocations, kg

rho_NTO = 1380; % Density of NTO @ 293K kg/m3
rho_MMH = 857; % Density of MMH @ 298K kg/m3

Prop_vol_ox = Prop_mass_ox./rho_NTO; % Propellant Oxidizer Volume Vector, m3
Prop_vol_fuel = Prop_mass_fuel./rho_MMH; % Propellant Fuel Volume Vector, m3

% Maximum propellant volume required m3
Prop_vol_ox_max = Prop_vol_ox(10);
Prop_vol_fuel_max = Prop_vol_fuel(10);

% Maximum Tank Volume Calculations (assuming 10% ullage) m3
Tank_vol_ox_max = Prop_vol_ox_max*(1+(10/100));
Tank_vol_fuel_max = Prop_vol_fuel_max*(1+(10/100));

% Ullage Volume Calculations m3
Ull_vol_ox = Tank_vol_ox_max - Prop_vol_ox_max;
Ull_vol_fuel = Tank_vol_fuel_max - Prop_vol_fuel_max;

%% Alternatively

% Requirement based on propellant mass (same parameters as above)

dV = 14.1979; % km/s
FM = 1000; % kg
c = 3.1889; % km/s

Prop_mass = FM*(exp(dV/c)-1); % kg

Prop_mass_flight_tot = Prop_mass*(1/0.85); % kg
Prop_mass_FC = 0.10*Prop_mass_flight_tot; % kg
Prop_mass_RP = 0.02*Prop_mass_flight_tot; % kg
Prop_mass_LU = 0.005*Prop_mass_flight_tot; % kg
Prop_mass_ONP = 0.01*Prop_mass_flight_tot; % kg
Prop_mass_ONO = 0.01*Prop_mass_flight_tot; % kg

Prop_mass_MM = (Prop_mass + Prop_mass_FC)*0.05; % kg
Prop_mass_cont = Prop_mass_flight_tot*0.05; % kg
Prop_mass_max = Prop_mass_flight_tot + Prop_mass_cont; % kg

Prop_mass_allo = [Prop_mass, Prop_mass_flight_tot, Prop_mass_FC, Prop_mass_RP, Prop_mass_LU, Prop_mass_ONP, Prop_mass_ONO, Prop_mass_MM, Prop_mass_cont, Prop_mass_max]; % kg

OFR = 2; % Oxidizer to fuel ratio
Prop_mass_ox = (OFR/(OFR+1)).*Prop_mass_allo; % Vector of oxidizer mass allocations, kg
Prop_mass_fuel = (1/(OFR+1)).*Prop_mass_allo; % Vector of fuel mass allocations, kg

rho_NTO = 1380; % Density of NTO @ 293K kg/m3
rho_MMH = 857; % Density of MMH @ 298K kg/m3

Prop_vol_ox = Prop_mass_ox./rho_NTO; % Propellant Oxidizer Volume Vector, m3
Prop_vol_fuel = Prop_mass_fuel./rho_MMH; % Propellant Fuel Volume Vector, m3

% Maximum propellant volume required m3
Prop_vol_ox_max = Prop_vol_ox(10);
Prop_vol_fuel_max = Prop_vol_fuel(10);

% Maximum Tank Volume Calculations (assuming 10% ullage) m3
Tank_vol_ox_max = Prop_vol_ox_max*(1+(10/100));
Tank_vol_fuel_max = Prop_vol_fuel_max*(1+(10/100));

% Ullage Volume Calculations m3
Ull_vol_ox = Tank_vol_ox_max - Prop_vol_ox_max;
Ull_vol_fuel = Tank_vol_fuel_max - Prop_vol_fuel_max;