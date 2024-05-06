# %%
import numpy as np
import math as math
import plotly.io as pio

from astropy import units as u
from astropy import time
from astropy.coordinates import solar_system_ephemeris

from poliastro import iod
from poliastro.bodies import Body, Earth, Jupiter, Sun
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.util import time_range
from poliastro.plotting import StaticOrbitPlotter

pio.renderers.default = "plotly_mimetype+notebook_connected"
solar_system_ephemeris.set("jpl")

# Initial data
date_launch = time.Time("2025-10-02 00:00", scale="utc").tdb
date_arrival = time.Time("2028-04-16 00:00", scale="utc").tdb

earth = Ephem.from_body(Earth, time_range(date_launch, end=date_arrival))
jupiter = Ephem.from_body(Jupiter, time_range(date_launch, end=date_arrival))
# Solve for departure and target orbits
ss_earth = Orbit.from_ephem(Sun, earth, date_launch)
ss_jupiter = Orbit.from_ephem(Sun, jupiter, date_arrival)

print(ss_earth.v/24/60/60)
print(ss_earth.r)
print(ss_jupiter.v/4/60/60)
print(ss_jupiter.r)

# Solve for the transfer maneuver
man_lambert = Maneuver.lambert(ss_earth, ss_jupiter)

# calculate dV magnitude
dv1_v = man_lambert.impulses[0][1]
dv2_v = man_lambert.impulses[1][1]
print (dv1_v)
print (dv2_v)
man_lambert.get_total_cost()

# Get the transfer and final orbits
ss_trans, ss_target = ss_earth.apply_maneuver(man_lambert, intermediate=True)

print(ss_trans.a)
print(ss_trans.ecc)

print(ss_target.a)
print(ss_target.ecc)

from poliastro.plotting import OrbitPlotter3D
plotter = OrbitPlotter3D()
plotter.set_attractor(Sun)

plotter.plot_ephem(earth, date_launch, label="Earth at launch position")
plotter.plot_ephem(jupiter, date_arrival, label="Jupiter at arrival position", color="purple")
plotter.plot_ephem(earth, date_arrival, label="Earth at arrival position")
plotter.plot_ephem(jupiter,date_launch,label="Jupiter at launch position", color="orange")
plotter.plot_trajectory(
    ss_trans.sample(max_anomaly=180 * u.deg),
    color="black",
    label="Transfer orbit",
)
plotter.set_view(45 * u.deg, 260 * u.deg, distance=3* u.km)

l1 = "Earth at launch (" + date_launch.strftime('%Y-%m-%d') + ")"
l2 = date_launch.strftime('%Y-%m-%d') + " (Jupiter at launch)"
l3 = date_arrival.strftime('%Y-%m-%d') + " (Earth at arrival)"
l4 = date_arrival.strftime('%Y-%m-%d') + " (Jupiter at arrival)"

plotter = StaticOrbitPlotter()
plotter.plot_body_orbit(
    Earth,
    date_launch,
    label="Earth at launch",
    color="C1",
    trail=True,
)
plotter.plot_ephem(jupiter,date_launch,label=l2, color="orange")
plotter.plot_ephem(earth, date_arrival, label=l3, color="blue")
plotter.plot_ephem(jupiter, date_arrival, label=l4, color="purple")
plotter.plot_trajectory(
    ss_trans.sample(max_anomaly=180 * u.deg),
    color="black",
    label="Transfer orbit",
)

# Jupiter to Europa transfer

from astropy.constants import Constant

my_mean_a_europa = Constant(
    "mean_a_europa",
    "My Europa's orbit mean semimajor axis around Jupiter",
    671100,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

my_GM_europa = Constant(
    "GM_europa",
    "My Europa gravitational constant",
    3202.71210,
    "km3 / (s2)",
    0.00050,
    "NASA Planetary Satellite Physical Parameters",
    system="si",
)

Europa = Body(
    parent=Jupiter,
    k=my_GM_europa,
    name="Europa",
    mass=4.8e22,
    R=2560,
    rotational_period=3.551*24*60*60
)

craft_earth = Orbit.circular(Earth,400 *u.km,epoch=date_launch)

# calculate excess hyperbolic speed first impulse
print(ss_trans.a)
print(ss_trans.ecc)
print((craft_earth.v + ss_trans.v))


# Calculate excess hyperbolic speed second impulse
craft_jupiter = Orbit.circular(Jupiter,1000000 *u.km,epoch=date_arrival)
print(ss_target.a)
print(ss_target.ecc)
print((craft_jupiter.v + ss_target.v))



