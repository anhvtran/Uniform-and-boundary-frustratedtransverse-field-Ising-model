# even system size 

from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
import json

def read_dict(dict_name):
    # for reading also binary mode is important
    with open(dict_name, 'rb') as fp:
        n_dict = json.load(fp)
        return n_dict

# Redefine:
read = read_dict("ed.json") # read the even section data 
epsilon_values = read['epsilon_values']
L_values = read['L_values']
h_values = read['h_values']
Egs_values = read['Egs_values'] # Actual global ground state 
First_excited_values = read['First_excited_values']
Second_excited_values = read['Second_excited_values']
Mag_values = read['Mag_values']

Egs = Egs_values[1]
First = First_excited_values[1]
Second = Second_excited_values[1]
color = 'cgbrmykcgbrmyk'

    ################# Plotting the Energy gaps ###############################
nfont = 18
plt.rcParams["font.size"] = nfont
plt.rcParams["font.family"] = "Arial"

### Gap between First Excited State and Ground State 

for i in  [3]:
    C = str(color[i])
    L_lab = 'L = '+ str(L_values[i])
    plt.plot(h_values,np.array(First[i])-np.array(Egs[i]), color=C, label = str(L_lab) )
    plt.legend(loc= "upper right", fontsize = nfont-4, framealpha=0.95)
    #plt.title("Energy gap between First Excited State and Ground State")
plt.xlabel("$\mathdefault{h/J}$")
plt.ylabel("$\mathdefault{(E_1 - E_0)/J}$")
plt.savefig('Gap1_10.png',dpi=300,bbox_inches = 'tight')
plt.show()

#   Plotting the magnetization varies on external field h 
line_values = ['solid','dotted','dashed']
for j in range(0,len(epsilon_values)):
    line = str(line_values[j])
    M = Mag_values[j]
    for i in  range(0,len(L_values)):
        C = str(color[i])
        L_lab = 'L = '+ str(L_values[i])
        plt.plot(h_values,M[i], color=C,ls = line, label = str(L_lab) if j == 0 else "" )
        plt.legend(loc= "upper right", fontsize=nfont-4)
        #plt.title("GS Magnetization vs External field")
    plt.xlabel("$\mathdefault{h/J}$")
    plt.ylabel("$\mathdefault{M}$")
plt.text(0.7, 0.5, "$\leftarrow$with negative $\mathdefault{\epsilon}$", fontsize=nfont-5)
plt.text(0.7, -0.5, "$\leftarrow$with positive $\mathdefault{\epsilon}$", fontsize=nfont-5)
plt.text(-0.075, -0.3, "without $\mathdefault{\epsilon}$", fontsize=nfont-5)
plt.arrow(0.07, -0.15, 0, 0.08, head_width=0.02, head_length=0.03, color='black')
plt.savefig('GS_Magnetization.png',dpi=300,bbox_inches = 'tight')
plt.show()