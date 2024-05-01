function [deltav1,deltav2,deltavtot,e, TOF] = hohmann(mu, periapsis, apoapsis)
% UNTITLED5 Summary of this function goes here
%   Hohmann transfer between two circular orbits
init

    % Hohmann transfers

    rpt = periapsis;            % periapsis transfer ellipse
    rat = apoapsis;              % apoapsis transfer  ellipse
    
    a_t = 0.5 * (rpt + rat);         % semi-major axis 
    v1i = v1_init(mu,rpt) ;          % initial velocity
    
    v2f = v2_final(mu,rat);        % final velocity    
    e = EccfromRpRa(rpt,rat);        % eccentricity of transfer orbit
    
    Et = EfromA(mu,a_t);           % Energy of transfer orbit
    vatransfer = VfromRA(mu,rat,a_t);  % velocity after delta v at destination
    vptransfer = VfromRA(mu,rpt,a_t);  % velocity after delta v at origin
    
    deltav1 = vptransfer - v1i;      % delta v origin
    deltav2 = v2f - vatransfer;      % delta v arrival
    deltavtot = abs(deltav1)+abs(deltav2) ;    % total delta v transfer
    TOF = 0.5 * PeriodfromA(mu,a_t);   % Time of flight transfer 
    TOFday = TOF/Day;
    TOFyear = TOF/Year;
    
    mmorigin = MMfromMua(mu,rpt);
    mmsc = MMfromMua(mu,a_t);      % mean motion of transfer 
    mmdestination = MMfromMua(mu,rat); % mean motion of destination
    synod1 = synodFromMm(mmsc,mmdestination);        % synodic period
    synodcheck = synodFromMua(mu,a_t,rat);     % verify synod
        
    synodyear = synod1/Year;
    
    phasei = pi -( mmdestination * TOF);  % phase angle for intercept 
    phasedegi = rad2deg(phasei);
    phasef = pi -( mmorigin * TOF);  % phase angle for intercept 
    phasedegf = rad2deg(phasef);
   
    
end

    


