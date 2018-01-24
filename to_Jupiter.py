
###############################################################################################################################
    
import numpy as np
from astropy import units as u
from astropy import time
from poliastro.bodies import Sun
from poliastro import iod

import pandas as p
import math
from poliastro.bodies import Earth, Mars, Jupiter
from poliastro.twobody import Orbit

import matplotlib.pyplot as plt
import warnings
from matplotlib.dates import DateFormatter
import datetime

###############################################################################################################################

# Wyznaczenie optymalnej daty tranzytu

def optimal_transit(date, transit_min, transit_max, planet1, planet2, vs0, step):

    
    date_arrival = date + transit_min   # minimalna data wykonania tranzytu
    date_max = date + transit_max       # maksymalna data wykonania, zakonczenie petli
    date_arrival_final = date_arrival

    vs_temp = 0 * u.km / u.s
    dv_final = 0 * u.km / u.s
    step_first = True

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

        dv_vector = vs1 - (vs0 + (v_1 / (24*3600) * u.day / u.s))           # zmiana predkosci niezbedna do udanego wykonania manewru
        dv = np.linalg.norm(dv_vector/10) * u.km / u.s                      # modul wektora zmiany predkosci

        if step_first:                                                        # zapis wynikow z pierwszego kroku
            dv_final = dv
            vs_temp = vs2

            step_first = False
        else:
            if dv < dv_final:                                               # sprawdzenie czy kolejny krok jest bardziej korzystna
                dv_final = dv
                date_arrival_final = date_arrival
                vs_temp = vs2

        date_arrival += step * u.day

    return dv_final, date_arrival_final, vs_temp                               # funkcja zwraca niezbedny przyrost predkosci, date przybycia 

###############################################################################################################################

# przypisanie minimalnych i maksymalnych dlugosci trwania manewru ziemia-mars oraz mars-jowisz

transit_minEM = 50 * u.day
transit_maxEM = 400 * u.day
transit_minMJ = 300 * u.day
transit_maxMJ = 1000 * u.day

def optimal_date(H, date0, date1, m, Isp, step):

    delta_v = 0 * u.km / u.s
    v_out = 0 * u.km / u.s
    m_prop = 0 * u.kg
    date_in = date0
    date_out = date0

    step_one0 = True
    x1 =[]
    x2 =[]
    x3 =[]
    x4 =[]
    
    # petla szukajaca daty startu z najmniejszym kosztem
    while date0 < date1:
        epoch0 = date0.jyear                            # przypisanie poczatkowej daty
        ss0 = Orbit.circular(Earth, H, epoch=epoch0)    # poczatkowa kolowa orbita dookola Ziemii
        vsE = ss0.rv()[1]                               # wektor predkosci na orbicie poczatkowej

        # optymalizacja kolejnych manewrow pod wzgledem daty zakonczenia tranzytu
        
        dvEM, date_arrivalM, vsM = optimal_transit(date0, transit_minEM, transit_maxEM, Earth, Mars, vsE, step)

        dvMJ, date_arrivalJ, vsJ = optimal_transit(date_arrivalM, transit_minMJ, transit_maxMJ, Mars, Jupiter, vsM, step)

        # koszty paliwa:
        m_pMJ = m * (math.exp((dvMJ) / Isp) - 1)
        m_pEM = (m + m_pMJ) * (math.exp((dvEM) / Isp) - 1)

        dv_total = (dvEM) + (dvMJ)                      # calkowita zmiana predkosci potrzebna do wykonania manewru
        m_p = m_pEM + m_pMJ                             # masa paliwa potrzebna do wykonania manewrow
        
        x1.append(date0.iso[0:10])
        x2.append(int((date_arrivalJ - date0).jd))
        x3.append(float(dv_total / u.km * u.s))
        x4.append(int(m_p / u.kg))
        #print(x1)
        
        #x={'1': date0.iso[0:10],'2': int((date_arrivalU - date0).jd),'3': float(dv_total / u.km * u.s),'4': int(m_p / u.kg)}
        #lista = pd.DataFrame(,index=[i])
        
        
        # wyswietlenie wynikow dla kolejnych dni
        print(date0.iso[0:10], ', %i days, %.3f km/s, %i kg' % (int((date_arrivalJ - date0).jd),float(dv_total / u.km * u.s),int(m_p / u.kg)))
 
        # zapisanie wynikow z pierwszego kroku
        if step_one0:
            delta_v = dv_total
            m_prop = m_p
            v_out = vsJ
            date_out = date_arrivalJ

            step_one0 = False
            # sprawdzenie czy nowa konfiguracja jest bardziej korzystna
        else:
            if dv_total < delta_v:
                delta_v = dv_total
                m_prop = m_p
                v_out = vsJ
                date_in = date0
                date_out = date_arrivalJ

        date0 += step * u.day
    listt = {'x1': x1,'x2': x2, 'x3': x3, 'x4': x4}
    listt = p.DataFrame(listt)
        
    return delta_v, v_out, date_in, date_out, m_prop, listt


#######################################################################################################################

def orbit_check(date, v, m, Isp):
    
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
        
##########################################################################################################################
        
# funkcja mająca na celu wylaczenie wyswietlania warningow
def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
    
#######################################################################################################################
    
    # Paremetry wejsciowe
     
    date_0 = '2018-01-01 12:00'
    date_1 = '2020-01-01 12:00'
    date_0 = time.Time(date_0, scale='utc')
    date_1 = time.Time(date_1, scale='utc')
    date_0 = time.Time(date_0.jd, format='jd', scale='utc')
    date_1 = time.Time(date_1.jd, format='jd', scale='utc')
    m_ship = 500 * u.kg
    I_sp = 500 * u.km / u.s
    H = 300 * u.km
    
    # czesc obliczeniowa programu
    
    step = 10 # krok czasowy symulacji w dniach
            
    # Asysta grawitacyjna Marsa (_M)
    
    delta_v_M, v_out_M, date_in_M, date_out_M, m_prop_M, list_M  = optimal_date(H, date_0, date_1, m_ship, I_sp, step)
    print (delta_v_M, v_out_M, date_in_M.iso, date_out_M.iso, m_prop_M)

    figure1, axis1 = plt.subplots()
    figure2, axis2 = plt.subplots()
    
    date_list = list_M.get('x1')
    dv_list = list_M.get('x3')
    mp_list = list_M.get('x4')
 
    # Rysowanie wykresow
    
    # wykres potrzebnej do wykonania manewru zmiany prędkosci
    
    axis1.grid(True)
    axis1.plot([datetime.datetime.strptime(val, '%Y-%m-%d') for val in date_list], dv_list, color='blue')
    
    figure1.autofmt_xdate()
    myFmt = DateFormatter("%Y-%m-%d")
    axis1.xaxis.set_major_formatter(myFmt)
    
    axis1.set_title('dV (start_date)')
    axis1.set_xlabel('start_date')
    axis1.set_ylabel('dV [km/s]')
    figure1.savefig("graph_dv.png")
    
    # wykres potrzebnej do wykonania manewru zmiany masy paliwa
    
    axis2.grid(True)
    axis2.plot([datetime.datetime.strptime(val, '%Y-%m-%d') for val in date_list], mp_list, color='red')
    
    figure2.autofmt_xdate()
    myFmt = DateFormatter("%Y-%m-%d")
    axis2.xaxis.set_major_formatter(myFmt)
    
    axis2.set_title('m_fuel (')
    axis2.set_xlabel('start_date')
    axis2.set_ylabel('mass_fuel [kg]')
    figure2.savefig("graph_mfuel.png")
    
    
