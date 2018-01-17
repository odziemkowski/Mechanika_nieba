import matplotlib.pyplot as plt
from astropy import units as u
from astropy import time
import warnings
from matplotlib.dates import DateFormatter
import datetime
import start_optimal_v2 as so
import chck_fun as c_f


# funkcja mająca na celu wylaczenie wyswietlania warningow
def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
    
#######################################################################################################################
    
    # Paremetry wejsciowe
     
    date_0 = '2020-01-01 12:00'
    date_1 = '2030-01-01 12:00'
    date_0 = time.Time(date_0, scale='utc')
    date_1 = time.Time(date_1, scale='utc')
    date_0 = time.Time(date_0.jd, format='jd', scale='utc')
    date_1 = time.Time(date_1.jd, format='jd', scale='utc')
    m_ship = 1000 * u.kg
    I_sp = 450 * u.km / u.s
    H = 30000 * u.km
    
    # czesc obliczeniowa programu
    
    step = 10 # krok czasowy symulacji w dniach
            
    # Asysta grawitacyjna Marsa
    delta_v_M, v_out_M, date_in_M, date_out_M, m_prop_M, lista_M  = so.start_date_optimal(H, date_0, date_1, m_ship, I_sp, step)
    print (delta_v_M, v_out_M, date_in_M.iso, date_out_M.iso, m_prop_M)
    c_f.check_jupiter_orbit(date_out_M, v_out_M, m_ship, I_sp)
    
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    
    daty_lista2 = lista_M.get('1')
    dv_lista2 = lista_M.get('3')
    mp_lista2 = lista_M.get('4')
 
    #Rysowanie wykresow predkosci i paliwa
    
    ax3.grid(True)
    ax3.plot([datetime.datetime.strptime(val, '%Y-%m-%d') for val in daty_lista2], dv_lista2, color='blue')
    
    fig3.autofmt_xdate()
    myFmt = DateFormatter("%Y-%m-%d")
    ax3.xaxis.set_major_formatter(myFmt)
    
    ax3.set_title('Wykres potrzebnej zmiany prędkosci od daty startu')
    ax3.set_xlabel('Data startu misji')
    ax3.set_ylabel('Całkowita zmiana prędkosci [km/s]')
    fig3.savefig("obraz3.png")
    # wykres 4
    ax4.grid(True)
    ax4.plot([datetime.datetime.strptime(val, '%Y-%m-%d') for val in daty_lista2], mp_lista2, color='green')
    
    fig4.autofmt_xdate()
    myFmt = DateFormatter("%Y-%m-%d")
    ax4.xaxis.set_major_formatter(myFmt)
    
    ax4.set_title('Ilosc paliwa potrzebna do realizacji misji')
    ax4.set_xlabel('Data startu misji')
    ax4.set_ylabel('Masa paliwa [kg]')
    fig4.savefig("obraz4.png")
    
    
    