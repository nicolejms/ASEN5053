clear
init
% generate date time arrays
t1Earth = datetime(2024,6,20,0,0,0);
t2Earth = datetime(2025,6,20,0,0,0);
tEarth = t1Earth:caldays(1):t2Earth;

TOF = PeriodfromA(Gsun,(AofEarthHelio+AofJupiterHelio)*AU)/2

% determine phase angles
    mmdestination = MMfromMua(Gsun,AofJupiterHelio*AU); % mean motion of destination    
    phasei = pi -( mmdestination * TOF);  % phase angle for intercept 
    phasedegi = rad2deg(phasei);

% create structures to store data
R= [0 0 0 0 0 0 0];Phase = {};
k = 1;
for i = 1:length(tEarth)
    [r1,~] = planetEphemeris(juliandate(tEarth(i)),'Sun','Earth',430);
    [r2,~] = planetEphemeris(juliandate(tEarth(i)),'Sun','Jupiter',430);     
    R = [R; juliandate(tEarth(i)) r1 r2];

    % assume both are on same ecliptic plane... so just x,y coord
    phase = acosd(dot(r1(1:2),r2(1:2))/(norm(r1(1:2))*norm(r2(1:2))));
    negphase = acosd(dot(r2(1:2),r1(1:2))/(norm(r1(1:2))*norm(r2(1:2))));
    if (phase > floor(-phasedegi) && phase < floor(-phasedegi+1)) 
        Phase = [Phase; string(tEarth(i))  phase];
    end
    if (negphase > floor(-phasedegi) && negphase < floor(-phasedegi+1)) 
        Phase = [Phase; string(tEarth(i))  negphase];
    end
end

