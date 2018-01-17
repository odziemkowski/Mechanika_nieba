import math
import numpy as np
from astropy import units as u
from astropy import time
from poliastro.bodies import Sun, Jupiter
from poliastro.twobody import Orbit

def check_jupiter_orbit(date, v, m, Isp):
    # sprawdzenie czy orbita została osiągnieta
    date_iso = time.Time(str(date.iso), format='iso', scale='utc')
    r_out=Orbit.from_body_ephem(Jupiter, date_iso)      # polozenie Jowisza po asyscie
    r_out1, vp_out1 = r_out.rv()
    v_exit = v + (vp_out1 / (24 * 3600) * u.day / u.s)       # predkosc satelity po manewrach
    epoch_out = date.jyear

    ss_out = Orbit.from_vectors(Sun, r_out1, v_exit, epoch=epoch_out)        # wyjsciowe parametry orbity 

    print('Sprawdzanie osiągnięcia orbity ...')
    print()

    if ss_out.ecc >= 1:     # ekscentrycznosc orbity
        print('Predkosc jest okej')
    else:
        print('Predkosc jest za mała')
        print()
        print('Dostosuj się do minimalnej orbity wyjściowej ')

        # minimalna orbita wyjsciowa paraboliczna:
        ss_out_new = Orbit.parabolic(Sun, ss_out.p, ss_out.inc, ss_out.raan, ss_out.argp, ss_out.nu, epoch=epoch_out)

        v_out_new = ss_out_new.rv()[1] - v_exit     # obliczenia nowej predkosci
        dv_out_new = np.linalg.norm(v_out_new) * u.km / u.s
        m_p_new = m * (math.exp(dv_out_new / Isp) - 1)        # obliczenia brakujace masy paliwa

        # Odpowiedz:
        print('Potrzebna delta V: %.3f km/s' % float(dv_out_new / u.km * u.s))
        print('Potrzebna dod. masa paliwa: %i kg' % int(m_p_new / u.kg))
        