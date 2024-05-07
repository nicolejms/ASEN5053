% This code will generate a pork chop plot of potential transfer orbits between planets.  
% It requires the lambert.m file found here: https://github.com/rodyo/FEX-Lambert
% Inputs area planet ids, departure timespan and target timespan in seconds
% since J2000. It also uses various functions from "Orbital Mechanics for
% Engineering Students" by Howard Curtis
% It outputs a pork chop plot with total delta-v requirements for an
% orbital transfer from planet to planet


clear all 
close all
init
% 'Pork Chop' plotter

% constants
mu = 1.32712440018e11; % [km^3/s^2]
dayInSec = 24*60*60;

% INPUTS

planetID1 = 3;
planetID2 = 5;

t1min = (J0(2025, 9, 10) - J0(2000,1,1))*dayInSec;
t1max = (J0(2025, 11, 10) - J0(2000,1,1))*dayInSec;
t2min = (J0(2028, 3, 1) - J0(2000,1,1))*dayInSec;
t2max = (J0(2028, 11, 1) - J0(2000,1,1))*dayInSec;

steps = 1000 ; % increasing increases computation time exponentially
deltaVlimit = 50; %[km/s] 

% MAIN CODE

% initialises matrices
x = zeros(1,steps);
y = zeros(1,steps);
z = zeros(1,steps);

counter = 1;

for t1 = t1min:((abs(t1max-t1min))/(steps-1)):t1max
    
    R1 = position(planetID1,t1);
    coe = coe_from_position(planetID1,t1);
    [r1, vPlanet1] = sv_from_coe(coe,Gsun);
    
    for t2 =  t2min:((abs(t2max-t2min))/(steps-1)):t2max
        
        R2 = position(planetID2,t2);
        coe = coe_from_position(planetID2,t2);
        [r2, vPlanet2] = sv_from_coe(coe,Gsun);
        
        tf = abs(t2-t1)/(dayInSec);
        [V1, V2, extremal_distances, exitflag] = lambert(R1, R2, tf, 0, mu);
        
        if exitflag == -1 % if lambert fails
            deltaV =  deltaVlimit+1; % ensures it is not displayed
            x(counter) = t1;
            y(counter) = t2;
            z(counter) = deltaV;
            counter = counter+1;             
            
        else
            deltaV = (abs(norm(vPlanet1-V1)))+(abs(norm(vPlanet2-V2))); % [km/s]
            x(counter) = t1;
            y(counter) = t2;
            z(counter) = deltaV;
            counter = counter+1;
            
        end
        
    end
    
end
launch = t1min/dayInSec;
arrival = t2min/dayInSec;

% finds lowest value of z
[M,I] = min(z);
minV = M
launcht = x(I) % [s]
arrivet = y(I) % [s]
minFueltf = abs(arrivet-launcht); % [s]

% comb out impossible deltaV's
for i=1:steps^2
    if z(i)>=deltaVlimit
        z(i)=NaN;
    else
    end
end

xv = linspace(min(x), max(x), steps);
yv = linspace(min(y), max(y), steps);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);

% 2D colourmap plot
figure(1)
peaks = [14.6 15 15.5 16 18 20 30 40 50 60]

contour(X,Y,Z,peaks,ShowText="on")

xticks(t1min:14*dayInSec:t1max);
datesX = {datestr(launch+(0*14),'dd-mmm-yy'),datestr(launch+(1*14),'dd-mmm-yy'),datestr(launch+(2*14),'dd-mmm-yy'),datestr(launch+(3*14),'dd-mmm-yy'),datestr(launch+(4*14),'dd-mmm-yy'),datestr(launch+(5*14),'dd-mmm-yy'),datestr(launch+(6*14),'dd-mmm-yy'),datestr(launch+(14*14),'dd-mmm-yy'),datestr(launch+(8*14),'dd-mmm-yy'),datestr(launch+(9*14),'dd-mmm-yy'),datestr(launch+(10*14),'dd-mmm-yy'),datestr(launch+(11*14),'dd-mmm-yy'),datestr(launch+(12*14),'dd-mmm-yy'),datestr(launch+(13*14),'dd-mmm-yy'),datestr(launch+(14*14),'dd-mmm-yy'),datestr(launch+(15*14),'dd-mmm-yy'),datestr(launch+(16*14),'dd-mmm-yy')};
xticklabels(datesX)
xtickangle(45)

yticks(t2min:30*dayInSec:t2max);
datesY = {datestr(arrival+(0*30),'dd-mmm-yy'), datestr(739983+(1*30),'dd-mmm-yy'),datestr(arrival+(2*30),'dd-mmm-yy'),datestr(arrival+(3*30),'dd-mmm-yy'),datestr(arrival+(4*30),'dd-mmm-yy'),datestr(arrival+(5*30),'dd-mmm-yy'),datestr(arrival+(6*30),'dd-mmm-yy'),datestr(arrival+(7*30),'dd-mmm-yy'),datestr(arrival+(8*30),'dd-mmm-yy'),datestr(arrival+(9*30),'dd-mmm-yy'),datestr(arrival+(10*30),'dd-mmm-yy'),datestr(arrival+(11*30),'dd-mmm-yy'), datestr(arrival+(12*30),'dd-mmm-yy'),datestr(arrival+(13*30),'dd-mmm-yy'),datestr(arrival+(14*30),'dd-mmm-yy'),datestr(arrival+(15*30),'dd-mmm-yy'),datestr(arrival+(16*30),'dd-mmm-yy'),datestr(arrival+(17*30),'dd-mmm-yy'),datestr(arrival+(18*30),'dd-mmm-yy'),datestr(arrival+(19*30),'dd-mmm-yy')};
yticklabels(datesY)
ytickangle(0)

xlabel("Earth Departure Date")
ylabel("Jupiter Arrival Date")
title("Total \DeltaV of Transfer Between Earth and Jupiter")
colormap turbo;

hold on
c = colorbar;
c.Label.String = '\DeltaV  [km/s]';

% plots crosshairs at optimal point
y1=get(gca,'ylim');
x1=get(gca,'xlim');
plot([launcht launcht],y1,'color','r');
plot(x1,[arrivet arrivet],'color','r');

