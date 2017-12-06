'''Projekt zaliczeniowy - Mechanika Nieba
Misja na Europę - księżyc Jowisza. 
Wersja: Testowa
Autor: Andrzej Odziemkowski'''

import math
import numpy as np
from astropy import units as u
from astropy import time

from poliastro import ephem
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.plotting import plot


date_launch = time.Time('2011-11-26 15:02', scale='utc')
date_arrival = time.Time('2012-08-06 05:17', scale='utc')
tof = date_arrival - date_launch
r0 = ephem.get_body_ephem(ephem.EARTH, date_launch)
r = ephem.planet_ephem(ephem.MARS, date_arrival)
from poliastro import iod
(v0, v), = iod.lambert(Sun.k, r0, r, tof)
v0