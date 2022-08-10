from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
import json

def read_dict(dict_name):
    # for reading also binary mode is important
    with open(dict_name, 'rb') as fp:
        n_dict = json.load(fp)
        return n_dict

# Redefine1:
read1 = read_dict("full_evenL_upto16_h0.7.json") # read the even section data
hz_values = read1['hz_values']
J_link_values = read1['J_link_values']
epsilon_values = read1['epsilon_values']

Leven_values = read1['L_values']
es0even_values = read1['es0_values'] # Actual global ground state 
es1even_values = read1['es1_values']
es2even_values = read1['es2_values']
mageven_values = read1['mag_values']
correven_values = read1['corr_values']
# Redefine2:
read2 = read_dict("full_oddL_upto17_h0.7.json") # read the odd section data 
Lodd_values = read2['L_values']
es0odd_values = read2['es0_values'] 
es1odd_values = read2['es1_values']
es2odd_values = read2['es2_values']
magodd_values = read2['mag_values']
corrodd_values = read2['corr_values']

# Plotting
color = 'cgbrmykcgbrmyk' # color codes 

es0even_h = es0even_values[1]
es1even_h = es1even_values[1]
es2even_h = es2even_values[1]

es0odd_h = es0odd_values[1]
es1odd_h = es1odd_values[1]
es2odd_h = es2odd_values[1]


nfont = 18
plt.rcParams["font.size"] = nfont
plt.rcParams["font.family"] = "Arial"
# Magnetization for even system sizes
line_values = ['solid','dotted','dashed']
for j in range(0,len(epsilon_values)):
    line = str(line_values[j])
    Meven = mageven_values[j]
    for i in  range(1,len(Leven_values)):
        C = str(color[i])
        L_lab = 'L = '+ str(Leven_values[i])
        Meven[i] = [round(num,0) for num in Meven[i]]
        plt.plot(J_link_values,Meven[i],color=C,ls = line, marker = str(i), ms = 12, markevery = 10, label = str(L_lab) if j == 0 else "" )
        plt.legend(loc='upper right', bbox_to_anchor=(1,0.955), fontsize=nfont-4)
        #plt.title("GS Magnetization when h = %s" %hz) 
    plt.xlabel("$\mathdefault{J_{link}/J}$")
    plt.ylabel("$\mathdefault{M_{stag}}$")
#plt.savefig('StagMag_full_evenL_upto16_h0.7.png',dpi=300,bbox_inches = 'tight')
plt.show()

# Magnetization for odd system sizes
line_values = ['solid','dotted','dashed']
for j in range(0,len(epsilon_values)):
    line = str(line_values[j])
    Modd = magodd_values[j]
    count = 1
    for i in [1, 2, 4, 5]:
        C = str(color[i])
        L_lab = 'L = '+ str(Lodd_values[i])
        plt.plot(J_link_values,Modd[i], color=C,ls = line, label = str(L_lab) if j == 0 else "" )
        plt.legend(loc='upper right', bbox_to_anchor=(1,0.955), fontsize=nfont-4)
        #plt.title("GS Magnetization when h = %s" %hz) 
        count += 1
    plt.xlabel("$\mathdefault{J_{link}/J}$")
    plt.ylabel("$\mathdefault{M_{stag}}$")
#plt.savefig('stagmag_full_oddL_upto17_h0.7.png',dpi=300,bbox_inches = 'tight')
plt.show()

# 1st gap at a given hz
for i in  range(0,len(Lodd_values)):
    L_lab = 'L = '+ str(Lodd_values[i])
    plt.plot(J_link_values,np.array(es1odd_h[i])-np.array(es0odd_h[i]),ls= 'solid',label = str(L_lab) )
for i in  range(0,len(Leven_values)):
    L_lab = 'L = '+ str(Leven_values[i])
    plt.plot(J_link_values,np.array(es1even_h[i])-np.array(es0even_h[i]),ls='dashed', label = str(L_lab) )
#plt.legend(loc='upper right', bbox_to_anchor=(1,0.90), prop={'size': 6})
plt.legend(loc=4, fontsize=nfont-5, framealpha=0.95)
plt.xlabel("$\mathdefault{J_{link}/J}$")
plt.ylabel("$\mathdefault{(E_1 - E_0)/J}$")
#plt.savefig('gap1_full_upto16_h0.7.png',dpi=300,bbox_inches = 'tight')
plt.show()

# 2nd gap at a given hz 

for i in  range(0,len(Lodd_values)):
    L_lab = 'L = '+ str(Lodd_values[i])
    plt.plot(J_link_values,np.array(es2odd_h[i])-np.array(es0odd_h[i]),ls= 'solid',label = str(L_lab) )
for i in  range(0,len(Leven_values)):
    L_lab = 'L = '+ str(Leven_values[i])
    plt.plot(J_link_values,np.array(es2even_h[i])-np.array(es0even_h[i]),ls='dashed', label = str(L_lab) )
#plt.legend(loc='upper right', bbox_to_anchor=(1,0.685), prop={'size': 6})
plt.legend(loc=4, fontsize=nfont-5, framealpha=0.95)
plt.xlabel("$\mathdefault{J_{link}/J}$")
plt.ylabel("$\mathdefault{(E_2 - E_0)/J}$")
#plt.savefig('gap2_full_upto16_h0.7.png',dpi=300,bbox_inches = 'tight')
plt.show()

# Plot the correlation relation 
correven_h = correven_values[1]
corrodd_h = corrodd_values[1]
for i in  range(0,len(Lodd_values)):
    L_lab = 'L = '+ str(Lodd_values[i])
    plt.plot(J_link_values,corrodd_h[i], ls = 'solid',label = str(L_lab) )
for i in  range(1,len(Leven_values)):
    L_lab = 'L = '+ str(Leven_values[i])
    plt.plot(J_link_values,correven_h[i],marker = str(i), ms = 12, markevery = 10, ls = 'dashed', label = str(L_lab) )
#plt.legend(loc='lower right',prop={'size': 6})
plt.legend(loc=4, fontsize=nfont-5, framealpha=0.95)
plt.xlabel("$\mathdefault{J_{link}/J}$")
plt.ylabel("$\mathdefault{C}$")
#plt.savefig('corr_full_evenL_upto16_h0.7.png',dpi=300,bbox_inches = 'tight')
plt.show()
