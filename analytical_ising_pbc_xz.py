import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import winsound

def Omega(k, h):
    return 2 * np.sqrt((np.cos(k)-h)**2 + np.sin(k)**2)

def Kplus1(L):
    n = np.arange(1, L//2)
    return 2*np.pi*n/L

def Kplus0(L):
    n = np.arange(1, L//2+1)
    return (2*n-1) * np.pi/L

# Execution
if __name__== '__main__':
    L = 2000

    start_time = datetime.now() # Start counting the time 
    hh = np.linspace(0, 2, 41)
    E00 = np.zeros(len(hh))
    E01 = np.zeros(len(hh))
    E10 = np.zeros(len(hh))
    E11 = np.zeros(len(hh))

    # Plotting 
    nfont = 18
    plt.rcParams["font.size"] = nfont
    plt.rcParams["font.family"] = "Arial"
    for ii in range(len(hh)):
        E00[ii] = -np.sum( Omega(Kplus0(L), hh[ii]) )
        E01[ii] = -2 - np.sum( Omega(Kplus1(L), hh[ii]) )
    plt.xlabel("$\mathdefault{h/J}$")
    plt.ylabel("$\mathdefault{(E_1 - E_0)/J}$")
    plt.plot(hh, E01-E00)
    plt.legend(["L = {first}".format(first=L)],loc= "upper right", fontsize = nfont-4, framealpha=0.95)
    plt.savefig('Gap1_2000_Analytic.png',dpi=300,bbox_inches = 'tight')
    plt.show()

    end_time = datetime.now()
    print("Calculation completed")
    print('Duration: {}'.format(end_time - start_time))
    winsound.Beep(500, 2500) # alarm when the code finished already

# or 
# import os 
# os.system('play -nq -t alsa synth {} sine {}'.format(5, 500))
# for Linux and Mac 
