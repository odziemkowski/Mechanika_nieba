'''Projekt zaliczeniowy - Mechanika Nieba
Misja na Europę - księżyc Jowisza. 
Wersja: Testowa
Autor: Andrzej Odziemkowski'''


import numpy as np
import matplotlib.pyplot as plt
plt.ion()

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris

from poliastro.bodies import Sun, Earth, Jupiter
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.iod import izzo
from poliastro.plotting import plot, OrbitPlotter
from poliastro.util import norm

solar_system_ephemeris.set("jpl")

# definiowanie Europy
#class _Europe(_Body):
#    parent = Jupiter
#    k = 3201.0
#    name = "Europe"
#    symbol = None
#    R = 1569.0



date_launch = Time("2011-08-05 16:25", scale='utc')     # przypisanie daty startu
date_arrival = Time("2016-07-05 03:18", scale='utc')    # przypisanie daty przybycia na orbitę Jowisza

# Initial state of the Earth
ss_e0 = Orbit.from_body_ephem(Earth, date_launch)       # przypisanie do ss_e0 orbity i pozycji Ziemi w dniu startu
r_e0, v_e0 = ss_e0.rv()                                 # wydobycie pozycji i predkosci dla Ziemi w dniu startu

op = OrbitPlotter()
op.plot(ss_e0)                                          # wykres orbity i pozycji Ziemi w dniu startu

# Wyswietlanie manewru i poczatkowej orbity

op = OrbitPlotter()
op.plot(ss_e0)

# rysowanie orbit
# ss_e0 to orbita Ziemii w chwili startu

op = OrbitPlotter()

op.plot(ss_e0, label="Earth launch state")

# And now, go to Jupiter!
ss_j = Orbit.from_body_ephem(Jupiter, date_arrival) # orbita i pozycja Jowisza w dni przybycia satelity do niego
r_j, v_j = ss_j.rv()                                # pozyskanie pozycji i prędkosc Jowisza w dniu przybucia satelity

# rozwiazanie problemu lamberta pomiedzy Ziemia i Jowiszem pomiedzy data przybycia i przelotu kolo Ziemii
(v_flypre, v_oip), = izzo.lambert(Sun.k, r_j, r_j, date_launch - date_arrival)    

# orbita wejscia na orbite jowisza
ss_oip = Orbit.from_vectors(Sun, r_j, v_oip, epoch=date_arrival)

fig, ax = plt.subplots(figsize=(9, 9))

op = OrbitPlotter(ax)

op.plot(ss_e0, label="Earth")
op.plot(ss_oip, label="Jupiter Orbit Insertion Phase")
op.plot(ss_j, label="Jupiter")