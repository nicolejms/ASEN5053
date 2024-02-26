init
g0 = 32.174;
It = 100;
H2gamma = 1.4;
H2isp = 284;
H2rho = 1.1;  % at 3500 psi/0 C
H2rho = 1.77;


m_prop = @(It,Isp) It/Isp/g0;
radius = @(vol) (vol/pi*0.75)^(1/3)

H_prop  = m_prop(It,H2isp)
H_vol = H_prop * H2rho
H_radius = radius(H_vol)


N2gamma = 1.4;
N2isp = 76;
N2rho = 24.7;
N2_prop = m_prop(It,N2isp)
N2_vol = N2_prop * N2rho
N2_radius = radius(N2_vol)

Argamma = 1.67;
Arisp = 57;
Arrho = 35.3;
Ar_prop = m_prop(It,Arisp)
Ar_vol = Ar_prop * Arrho
Ar_radius = radius(Ar_vol)
