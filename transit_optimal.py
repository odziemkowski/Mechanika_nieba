import numpy as np
from astropy import units as u
from astropy import time
from poliastro.bodies import Sun
from poliastro.twobody import Orbit
from poliastro import iod

# Wyznaczenie optymalnej daty zakonczenia tranzytu
def transit_optimal(date, transit_min, transit_max, planet1, planet2, vs0, step):

    
    date_arrival = date + transit_min   # minimalna data tranzytu
    date_max = date + transit_max       # maksymalna data, zakonczenie petli
    date_arrival_final = date_arrival

    vs22 = 0 * u.km / u.s
    dv_final = 0 * u.km / u.s
    step_one = True

    # petla idaca po datach z okreslonym krokiem
    while date_arrival < date_max:
        tof = date_arrival - date                                           # tof - time of flight
        date_iso = time.Time(str(date.iso), scale='utc')                    # data startu 
        date_arrival_iso = time.Time(str(date_arrival.iso), scale='utc')    # data przylotu

        r1= Orbit.from_body_ephem(planet1, date_iso)
        r2= Orbit.from_body_ephem(planet2, date_arrival_iso)
        r_1, v_1 = r1.rv()                                                  # pozycja i predkosc planety poczatkowej
        r_2, v_2 = r2.rv()                                                  # pozycja i predkosc planety koncowej
        (vs1, vs2), = iod.lambert(Sun.k, r_1, r_2, tof)                     # rozwiazanie zagadnienia lamberta

        dv_vector = vs1 - (vs0 + (v_1 / (24*3600) * u.day / u.s))           # niezbedna do udanego wykonania manewru zmiana predkosci
        dv = np.linalg.norm(dv_vector/10) * u.km / u.s                      # modul wektora zmiany predkosci

        if step_one:                                                        # zapis wynikow z pierwszego kroku
            dv_final = dv
            vs22 = vs2

            step_one = False
        else:
            if dv < dv_final:                                               # sprawdzenie czy kolejny krok jest bardziej korzystna
                dv_final = dv
                date_arrival_final = date_arrival
                vs22 = vs2

        date_arrival += step * u.day

    return dv_final, date_arrival_final, vs22                               # funkcja zwraca niezbedny przyrost predkosci, date przybycia 
