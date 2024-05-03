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

%maxtf = 250*dayInSec; % [s]

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
% pcolor(X,Y,Z) % effectivly a 2D surf
% shading interp

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

% datestr(737855) is opposite of datenum('2020,3,5')
% datenum('2000,0,0') = 730485

hold on
c = colorbar;
c.Label.String = '\DeltaV  [km/s]';

% plots crosshairs at optimal point
y1=get(gca,'ylim');
x1=get(gca,'xlim');
plot([launcht launcht],y1,'color','r');
plot(x1,[arrivet arrivet],'color','r');

% jan 1 2029
% dec 10 2031
%%%%%%%%%%%%%%%%%%
% other plots
% figure(4)
% surf(X,Y,Z);
% plot3(x,y,z)
% stem3(x,y,z)
%%%%%%%%%%%%%%%%%%%

% 
% % calculates difference between NASA prediction and the software's
% departureDays = launcht/dayInSec;
% NASAdeparture = (J0(2020, 7, 24.6) - J0(2000,1,1));
% differenceInDeparture = departureDays-NASAdeparture
% arrivetDays = arrivet/dayInSec;
% NASAarrive = (J0(2021, 2, 18) - J0(2000,1,1));
% differenceInArrival = arrivetDays-NASAarrive
% 
% % plots NASA launch window
% tNASAdepartureMin = (J0(2020, 7, 18) - J0(2000,1,1))*dayInSec;
% tNASAdepartureMax = (J0(2020, 8, 5) - J0(2000,1,1))*dayInSec;
% % plot([tNASAdepartureMin tNASAdepartureMax],[NASAarrive*dayInSec NASAarrive*dayInSec],'color','m');
% % hold off
% 
% %variation in deltaV for NASA window
% figure(2)
% R2 = position(planetID2,NASAarrive*dayInSec);
% coe = coe_from_position(planetID2,NASAarrive*dayInSec);
% [r2, vPlanet2] = sv_from_coe(coe);
% 
% hold on
% deltaVs = [];
% idc = 1;
% for tlaunch = tNASAdepartureMin-(24*60*60) : 12*60*60 : tNASAdepartureMax-(24*60*60)
% 
%     R1 = position(planetID1,tlaunch);
%     coe = coe_from_position(planetID1,tlaunch);
%     [r1, vPlanet1] = sv_from_coe(coe);
%     tf = abs(tlaunch-(NASAarrive*dayInSec))/(dayInSec);
%     [V1, V2, extremal_distances, exitflag] = lambert(R1, R2, tf, 0, mu);
%     deltaV = (abs(norm(vPlanet1-V1)))+(abs(norm(vPlanet2-V2))); % [km/s]
%     plot(tlaunch, deltaV,'.','color','b')
%     deltaVs(idc) = deltaV;
%     idc =idc+1;
%     % attempt at plotting my departure date
% %     y1=get(gca,'ylim');
% %     plot([launcht launcht],[6:7],'color','r');
% 
% end
% 
% [value,position] = min(deltaVs);
% mindVNASALaunchWindow = value;
% 
% % sets axes to dates
% xticks(tNASAdepartureMin-(24*60*60):24*60*60:tNASAdepartureMax-(24*60*60));
% datesX = {datestr(737990+(0*1),'dd-mmm-yy'),datestr(737990+(1*1),'dd-mmm-yy'),datestr(737990+(2*1),'dd-mmm-yy'),datestr(737990+(3*1),'dd-mmm-yy'),datestr(737990+(4*1),'dd-mmm-yy'),datestr(737990+(5*1),'dd-mmm-yy'),datestr(737990+(6*1),'dd-mmm-yy'),datestr(737990+(7*1),'dd-mmm-yy'),datestr(737990+(8*1),'dd-mmm-yy'),datestr(737990+(9*1),'dd-mmm-yy'),datestr(737990+(10*1),'dd-mmm-yy'),datestr(737990+(11*1),'dd-mmm-yy'),datestr(737990+(12*1),'dd-mmm-yy'),datestr(737990+(13*1),'dd-mmm-yy'),datestr(737990+(14*1),'dd-mmm-yy'),datestr(737990+(15*1),'dd-mmm-yy'),datestr(737990+(16*1),'dd-mmm-yy'),datestr(737990+(17*1),'dd-mmm-yy'),datestr(737990+(18*1),'dd-mmm-yy')};
% xticklabels(datesX)
% xtickangle(45)
% %title({'Variation in \DeltaV over NASAs proposed launch','window assuming an arrival date of 18-02-21'})
% xlabel('Departure Date')
% ylabel('\DeltaV  [km/s]')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % variation in deltaV for NASA window but my arrival date
% 
% % R2 = position(planetID2,arrivet);
% % coe = coe_from_position(planetID2,arrivet);
% % [r2, vPlanet2] = sv_from_coe(coe);
% 
% 
% % for tlaunch = tNASAdepartureMin-(24*60*60*4) : 60*60 : tNASAdepartureMax
% %     
% %     R1 = position(planetID1,tlaunch);
% %     coe = coe_from_position(planetID1,tlaunch);
% %     [r1, vPlanet1] = sv_from_coe(coe);
% %     tf = abs(tlaunch-(NASAarrive*dayInSec))/(dayInSec);
% %     [V1, V2, extremal_distances, exitflag] = lambert(R1, R2, tf, 0, mu);
% %     deltaV = (abs(norm(vPlanet1-V1)))+(abs(norm(vPlanet2-V2))); % [km/s]
% %     plot(tlaunch, deltaV,'.','color','m')
% %     
% % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% hold off
% 

 julian = 2451545.5
 launch = julian + launcht/Day
 arrival = julian + arrivet/Day
% 777527000
