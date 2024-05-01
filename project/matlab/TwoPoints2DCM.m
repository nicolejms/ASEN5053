function Rab = TwoPoints2DCM(ra, rb, va, vb)
% TWOPOINTS2DCM(ra, rb, va, vb) calculates the DCM, Rab, between the two
% reference frames, a & b, where vec_a = Rab * vec_b. The function requires
% the knowledge of the representations of two arbitrary points (vectors) 
% in both frames.
%
% INPUTS:
%   - ra: vector r represented in frame a
%   - rb: vector r represented in frame b
%   - va: vector v represented in frame a
%   - vb: vector v represented in frame b
% 
% OUTPUTS:
%   - Rab: An orthonaormal DCM that rotates any vector from frame b to
%          frame a
xac = normalize(ra, 'norm', 2);
xbc = normalize(rb, 'norm', 2);
zac = normalize(cross(ra, va), 'norm',2);
zbc = normalize(cross(rb, vb), 'norm',2);
yac = cross(zac, xac);
ybc = cross(zbc, xbc);
Rac = [xac; yac; zac];
Rbc = [xbc; ybc; zbc];
Rab = Rac*Rbc';
end