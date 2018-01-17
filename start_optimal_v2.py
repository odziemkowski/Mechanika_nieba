import pandas as pd
import math
from astropy import units as u
from poliastro.bodies import Earth, Mars, Jupiter
from poliastro.twobody import Orbit


import transit_optimal as t_o

# przypisanie minimalnych i maksymalnych dlugosci trwania manewru ziemia-mars oraz mars-jowisz

transit_minEM = 50 * u.day
transit_maxEM = 100 * u.day
transit_minMJ = 400 * u.day
transit_maxMJ = 700 * u.day

def start_date_optimal(H, date0, date1, m, Isp, step):

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
        
        dvEM, date_arrivalM, vsM = t_o.transit_optimal(date0, transit_minEM, transit_maxEM, Earth, Mars, vsE, step)

        dvMJ, date_arrivalJ, vsJ = t_o.transit_optimal(date_arrivalM, transit_minMJ, transit_maxMJ, Mars, Jupiter, vsM, step)

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
            # sprawdzenei czy nowa konfiguracja jest bardziej korzystna
        else:
            if dv_total < delta_v:
                delta_v = dv_total
                m_prop = m_p
                v_out = vsJ
                date_in = date0
                date_out = date_arrivalJ

        date0 += step * u.day
    lista = {'1': x1,'2': x2, '3': x3, '4': x4}
    lista = pd.DataFrame(lista)
        
    return delta_v, v_out, date_in, date_out, m_prop, lista
