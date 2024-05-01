% planet_plot.m      January 30, 2022

% graphics display of the orbital characteristics of the planets

% JPL DE440s ephemeris and MICE routines

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global iephem ephname km

% astronomical unit (kilometers)

aunit = 149597870.691;

% sun gravitational constant (kilometers^3/second^2)

smu = 132712440040.9446;

% angular conversion factors

dtr = pi / 180.0;

rtd = 180.0 / pi;

% define name of ephemeris binary file

ephname = 'de440s.bsp';

% initialize JPL ephemeris

iephem = 1;

% request output in kilometers and kilometers/second

km = 1;

% begin simulation

clc; home;

fprintf('\nprogram planet_plot - orbital characteristics of the planets\n');

while(1)
    
   fprintf('\nplanet menu\n');

   fprintf('\n  <1> Mercury\n');
   fprintf('\n  <2> Venus\n');
   fprintf('\n  <3> Earth\n');
   fprintf('\n  <4> Mars\n');
   fprintf('\n  <5> Jupiter\n');
   fprintf('\n  <6> Saturn\n');
   fprintf('\n  <7> Uranus\n');
   fprintf('\n  <8> Neptune\n');
   fprintf('\n  <9> Pluto');

   fprintf('\n\nplease select the planet\n');
          
   iplanet = input('? ');
            
   if (iplanet >= 1 && iplanet <= 9)
       
      break;
      
   end
   
end

% set name of target body

switch iplanet
    
    case (1)
        
        targ = 'Mercury';
        
    case (2)
        
        targ = 'Venus';
        
    case (3)
        
        targ = 'Earth';
        
    case (4)
        
        targ = 'Mars';
        
    case(5)
        
        targ = 'Jupiter';
        
    case (6)
        
        targ = 'Saturn';
        
    case (7)
        
        targ = 'Uranus';
        
    case (8)
        
        targ = 'Neptune';
        
    case (9)
        
        targ = 'Pluto';
        
end

% request initial calendar date

[month, day, year] = getdate;

while (1)
    
    fprintf('\nplease input the simulation period (days)\n');
    
    ndays = input('? ');
    
    if (ndays > 0)
        
        break;
        
    end
    
end

while (1)
    
    fprintf('\nplease input the graphics step size (days)\n');
    
    dtstep = input('? ');
    
    if (dtstep > 0)
        
        break;
        
    end
    
end

while (1)
    
    fprintf('\nplease select the orbital element to plot\n');
    
    fprintf('\n  <1> semimajor axis');
    
    fprintf('\n\n  <2> eccentricity');
    
    fprintf('\n\n  <3> orbital inclination');
    
    fprintf('\n\n  <4> argument of perihelion');
    
    fprintf('\n\n  <5> longitude of the ascending node');
    
    fprintf('\n\n  <6> true anomaly');
    
    fprintf('\n\n  <7> argument of latitude');
    
    fprintf('\n\n  <8> aphelion radius');
    
    fprintf('\n\n  <9> perihelion radius');
    
    fprintf('\n\n  <10> heliocentric distance');
    
    fprintf('\n\n  <11> heliocentric speed\n\n');
       
    oetype = input('? ');
    
    if (oetype >= 1 && oetype <= 11)
        
        break;
        
    end
    
end

fprintf('\n\n');

% compute initial tdb julian date

jdatei = julian(month, day, year);

t = -dtstep;

npts = ndays / dtstep;

% create and open ascii csv state vector data file

filename = horzcat(targ, '_state.csv');

fid1 = fopen(filename, 'w');

% create and open ascii csv orbital elements data file

filename = horzcat(targ, '_orbit.csv');

fid2 = fopen(filename, 'w');

% allocate x-y data points

x = zeros(npts, 1);

y = zeros(npts, 1);

for i = 1:1:npts
    
    t = t + dtstep;
    
    % current TDB julian date
    
    tjd = jdatei + t;
    
    % compute EME2000 state vector (kilometers and kilometers/second)
    
    rrd = jpleph_mice(tjd, iplanet, 11);
    
    % extract position and velocity vectors
    
    r1(1) = rrd(1);
    r1(2) = rrd(2);
    r1(3) = rrd(3);
    
    v1(1) = rrd(4);
    v1(2) = rrd(5);
    v1(3) = rrd(6);
    
    oev1 = eci2orb1(smu, r1, v1);
    
    % write TDB julian date and state vector to data file
    
    fprintf(fid1, '%+16.14e, %+16.14e, %+16.14e, %+16.14e, %+16.14e, %+16.14e, %+16.14e\n', ...
        tjd, r1(1), r1(2), r1(3), v1(1), v1(2), v1(3));
    
    % write TDB julian date and orbital elements to data file
    
    fprintf(fid2, '%+16.14e, %+16.14e, %+16.14e, %+16.14e, %+16.14e, %+16.14e, %+16.14e\n', ...
        tjd, oev1(1), oev1(2), rtd * oev1(3), rtd * oev1(4), rtd * oev1(5), rtd * oev1(6));
        
    % determine serial date number
    
    [month, day, year] = gdate(tjd);
    
    % load serial date number
    
    x(i) = datenum(year, month, day);
    
    switch oetype
        
        case 1
            
            % semimajor axis (AU)
            
            y(i) = oev1(1) / aunit;
            
        case 2
            
            % orbital eccentricity
            
            y(i) = oev1(2);
            
        case 3
            
            % orbital inclination (degrees)
            
            y(i) = oev1(3) * rtd;
            
        case 4
            
            % argument of perihelion (degrees)
            
            y(i) = oev1(4) * rtd;
            
        case 5
            
            % longitude of ascending node (degrees)
            
            y(i) = oev1(5) * rtd;
            
        case 6
            
            % true anomaly (degrees)
            
            y(i) = oev1(6) * rtd;
            
        case 7
            
            % argument of latitude (degrees)
            
            y(i) = mod(oev1(5) + oev1(6), 2.0 * pi) * rtd;
            
        case 8
            
            % aphelion radius (AU)
            
            y(i) = (oev1(1) * (1.0 + oev1(2))) / aunit;
            
        case 9
            
            % perihelion radius (AU)
            
            y(i) = (oev1(1) * (1.0 - oev1(2))) / aunit;
            
        case 10
            
            % heliocentric distance (AU)
            
            y(i) = norm(r1) / aunit;
            
        case 11
            
            % heliocentric speed (kilometers/second)
            
            y(i) = norm(v1);
              
    end
    
end

% close ascii csv data files

fclose(fid1);

fclose(fid2);

% create and label plot

figure(1);

plot(x, y, '-r', 'LineWidth', 1.5);

% label x axis with calendar dates

datetick('x', 'mm/dd/yyyy');

% define main plot and y-axis labels

switch oetype
    
    case 1
        
        main_title = horzcat('Semimajor Axis of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('semimajor axis (AU)', 'FontSize', 12);
        
    case 2
        
        main_title = horzcat('Orbital Eccentricity of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('orbital eccentricity', 'FontSize', 12);
        
    case 3
        
        main_title = horzcat('Orbital Inclination of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('orbital inclination (degrees)', 'FontSize', 12);
        
    case 4
        
        main_title = horzcat('Argument of Perihelion of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('argument of perihelion (degrees)', 'FontSize', 12);
        
    case 5
        
        main_title = horzcat('Longitude of Ascending node of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('longitude of ascending node (degrees)', 'FontSize', 12);
        
    case 6
        
        main_title = horzcat('True Anomaly of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('true anomaly (degrees)', 'FontSize', 12');
        
    case 7
        
        main_title = horzcat('Argument of Latitude of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('argument of latitude (degrees)', 'FontSize', 12');
        
    case 8
        
        main_title = horzcat('Aphelion Radius of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('aphelion radius (AU)', 'FontSize', 12');
        
    case 9
        
        main_title = horzcat('Perihelion Radius of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('perihelion radius (AU)', 'FontSize', 12');
        
    case 10
        
        main_title = horzcat('Heliocentric Distance of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('heliocentric distance (AU)', 'FontSize', 12');
        
    case 11
        
        main_title = horzcat('Heliocentric Speed of ', targ);
        
        title (main_title, 'FontSize', 16);
        
        ylabel('heliocentric speed (kilometers/second)', 'FontSize', 12');
                
end

xlabel('calendar date', 'FontSize', 12);

grid;

zoom on;

% create color eps graphics image with tiff preview

filename = horzcat(targ, '_plot.eps');

print('-depsc', filename);

% close all files

fclose('all');

cspice_kclear;

