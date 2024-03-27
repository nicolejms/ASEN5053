%% Problem 2
clear; clc;
close all;
init;
mu = Gsun;
path = "short";

% import data files 
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ["\t", ","];

% Specify column names and types
opts.VariableNames = ["JDTDB", "CalendarDateTDB", "X", "Y", "Z", "VX", "VY", "VZ", "VarName9"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["CalendarDateTDB", "VarName9"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["CalendarDateTDB", "VarName9"], "EmptyFieldRule", "auto");
%filename: "/Users/nicolejames/src/asen5050/asen5050/hwk/hwk6/HW6_Ephem_Earth.txt"
% Import the data

HW6EphemEarth = readtable("earth_eph.txt", opts);
HW6EphemMars = readtable("jupiter_eph_week.txt", opts);


R1 = [HW6EphemEarth.X, HW6EphemEarth.Y, HW6EphemEarth.Z];
V1 = [HW6EphemEarth.VX, HW6EphemEarth.VY, HW6EphemEarth.VZ];
R2 = [HW6EphemMars.X, HW6EphemMars.Y, HW6EphemMars.Z];
V2 = [HW6EphemMars.VX, HW6EphemMars.VY, HW6EphemMars.VZ];

% create structures to store data
dvE = []; dvM = []; TOFf = [];Epochs = [];
k = 1;
for i = 1:length(R1)
    r1 = R1(i,1:3);
    v1 = V1(i,1:3);
    for j = 1:length(R2)
        TOF = (HW6EphemMars.JDTDB(j) - HW6EphemEarth.JDTDB(i))*Day;
        if TOF<0
            disp("invalid TOF");
            break;    
        end
        r2 = R2(j,1:3);
        v2 = V2(j,1:3);
        [vinfe, vinfm] = lambert2(mu,r1,r2,v1,v2,TOF,path);

        VinfE(i,j) = norm(vinfe);
        VinfM(i,j) = norm(vinfm);
    end
end


X = HW6EphemEarth.JDTDB(1:end) - HW6EphemEarth.JDTDB(1);   
Y = HW6EphemMars.JDTDB(1:end) - HW6EphemMars.JDTDB(1); 

% X = datetime(HW6EphemEarth.JDTDB,'convertfrom','juliandate')
% Y = datetime(HW6EphemMars.JDTDB,'convertfrom','juliandate')
figure(1) 
% peaks = VinfE(:); 
% peaks = fix(peaks);

contour(X(1:536),Y(1:366),VinfE.',peaks,'ShowText','on')
%peaks = [3.9 4 4.1 4.2 4.3 4.5 4.7 5 5.5 6 7 8 13 20 30 45 60 61 62 63];

xlabel("Earth Departure (Days past 6/20/2005)")
ylabel("Mars Arrival (Days past 12/1/2005)")
title("V-infinity of Earth Departures")
colormap winter;
legend("V_\infty (km/s)", 'Location', 'northoutside','FontSize', 12)
axis equal;
set(gca,'XMinorTick','on','YMinorTick','on')


figure(2) 
% peaks = VinfM(:);
% peaks = fix(peaks);
%peaks = [2 2.4 2.5 2.6 2.7 2.8 2.9 3 3.5 3.8 4 4.5 5 6 7 8 20 30 35 41 42  43 44];
contour(X(1:748),Y,VinfM.',peaks,'ShowText','on')
xlabel("Earth Departure (Days past 6/20/2005)")
ylabel("Mars Arrival (Days past 12/1/2005)")
title("V-infinity of Mars Arrivals")
colormap turbo;
legend("V_\infty (km/s)", 'Location', 'northoutside','FontSize', 12)
axis equal;
set(gca,'XMinorTick','on','YMinorTick','on')



