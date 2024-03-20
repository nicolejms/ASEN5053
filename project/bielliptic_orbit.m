function deltavtot = bielliptic_orbit(mu, initial_orbit, intermediate_orbit, final_orbit)
% UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%%  bi elliptic transfer
init
    % transfer 1
    % rb = bielliptic_radius;
    % r1i = AofEarthHelio*AU;
    v1i = v1_init(mu,initial_orbit) % initial velocity of satellite
    v3f = v2_final(mu,final_orbit)
     
    at1 = 0.5*(initial_orbit + intermediate_orbit)       % semi major axis transfer 1
    et1 = (intermediate_orbit - initial_orbit)/(intermediate_orbit + initial_orbit) % eccentricity transfer 1
    
    % transfer 2
    at2 = 0.5*(intermediate_orbit + final_orbit)          % semi major axis transfer 2
    et2 = EccfromRpRa(intermediate_orbit, final_orbit)   % eccentricity transfer 2
    
    vpt1 = VfromRA(mu,initial_orbit,at1)    % velocity at periapsis transfer 1
    vat1 = VfromRA(mu,intermediate_orbit,at1)     % velocity at apoapsis transfer 1
    
    vpt2 = VfromRA(mu,final_orbit,at2)    % velocity at periapsis transfer 2
    vat2 = VfromRA(mu,intermediate_orbit,at2)     % velocity at apoapsis transfer 2
    deltav1t1 = vpt1 - v1i
    deltav2t1 = vat2 - vat1
    deltav3t2 = v3f - vpt2

deltavb = abs(deltav1t1) + abs(deltav2t1) + abs(deltav3t2)

TOFbi = pi * (sqrt(at1^3/Gsun) + sqrt(at2^3)/Gsun)
TOFbiYear = TOFbi/Year


% part 
deltavtot = deltavb
mratio = Mfinal(1,deltavtot,212)


end