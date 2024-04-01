  %% Problem 2
clear; clc;
close all;
init;
mu = Gsun;
path = "short";
% generate date time arrays
t1Earth = datetime(2024,6,1,0,0,0);
t2Earth = datetime(2028,6,1,0,0,0);
tEarth = t1Earth:caldays(2):t2Earth;

t1Jupiter = datetime(2026,12,1,0,0,0);
t2Jupiter = datetime(2031,12,1,0,0,0);
tJupiter = t1Jupiter:caldays(4):t2Jupiter;

% R1 = [HW6EphemEarth.X, HW6EphemEarth.Y, HW6EphemEarth.Z];
% V1 = [HW6EphemEarth.VX, HW6EphemEarth.VY, HW6EphemEarth.VZ];
% R2 = [HW6EphemMars.X, HW6EphemMars.Y, HW6EphemMars.Z];
% V2 = [HW6EphemMars.VX, HW6EphemMars.VY, HW6EphemMars.VZ];

% create structures to store data
dvE = []; dvM = []; TOFf = [];Epochs = [];R1 = [];R2 = [];V1 = []; V2 = [];
k = 1;
for i = 1:length(tEarth)
    [r1,v1] = planetEphemeris(juliandate(tEarth(i)),'Sun','Earth',430);
    R1 = [R1;r1];
    V1 = [V1;v1];
    for j = 1:length(tJupiter)
        TOF = juliandate(tJupiter(j)) - juliandate(tEarth(i));
        if TOF<0
            disp("invalid TOF");
            break;    
        end
        [r2,v2] = planetEphemeris(juliandate(tJupiter(i)),'Sun','Jupiter',430);
        R2 = [R2;r2];
        V2 = [V2;v2];
        [vinfe, vinfm] = lambert2(mu,r1,r2,v1,v2,TOF,path);

        VinfE(i,j) = norm(vinfe);
        VinfM(i,j) = norm(vinfm);
    end
end


X = R1.JDTDB(1:end) - HW6EphemEarth.JDTDB(1);   
Y = HW6EphemMars.JDTDB(1:end) - HW6EphemMars.JDTDB(1); 


figure(1) 
% peaks = VinfE(:); 
% peaks = fix(peaks);

contour(X(1:748),Y(1:422),VinfE.','ShowText','on')
%peaks = [3.9 4 4.1 4.2 4.3 4.5 4.7 5 5.5 6 7 8 13 20 30 45 60 61 62 63];

xlabel("Earth Departure (Days past 3/14/2024)")
ylabel("Jupiter Arrival (Days past 3/26/2027)")
title("V-infinity of Earth Departures")
colormap winter;
legend("V_\infty (km/s)", 'Location', 'northoutside','FontSize', 12)
% axis equal;
set(gca,'XMinorTick','on','YMinorTick','on')


figure(2) 
% peaks = VinfM(:);
% peaks = fix(peaks);
%peaks = [2 2.4 2.5 2.6 2.7 2.8 2.9 3 3.5 3.8 4 4.5 5 6 7 8 20 30 35 41 42  43 44];
contour(X(1:748),Y(1:422),VinfM.','ShowText','on')
xlabel("Earth Departure (Days past 3/14/2024)")
ylabel("Jupiter Arrival (Days past 3/26/2027)")
title("V-infinity of Jupiter Arrival")
colormap turbo;
legend("V_\infty (km/s)", 'Location', 'northoutside','FontSize', 12)
%axis equal;
set(gca,'XMinorTick','on','YMinorTick','on')



