# %%
from astropy import units as u
from astropy import time

import numpy as np

from poliastro import iod
from poliastro.bodies import Earth, Jupiter, Sun
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.util import time_range

# %%
# More info: https://plotly.com/python/renderers/
import plotly.io as pio

pio.renderers.default = "plotly_mimetype+notebook_connected"

# %%
from astropy.coordinates import solar_system_ephemeris

solar_system_ephemeris.set("jpl")

# %%
# Initial data
date_launch = time.Time("2024-08-17 15:02", scale="utc").tdb
date_arrival = time.Time("2027-05-11 05:17", scale="utc").tdb

# %%
earth = Ephem.from_body(Earth, time_range(date_launch, end=date_arrival))
jupiter = Ephem.from_body(Jupiter, time_range(date_launch, end=date_arrival))
# Solve for departure and target orbits
ss_earth = Orbit.from_ephem(Sun, earth, date_launch)
ss_jupiter = Orbit.from_ephem(Sun, jupiter, date_arrival)

# %%
# Solve for the transfer maneuver
man_lambert = Maneuver.lambert(ss_earth, ss_jupiter)

# Get the transfer and final orbits
ss_trans, ss_target = ss_earth.apply_maneuver(man_lambert, intermediate=True)

# %%
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


