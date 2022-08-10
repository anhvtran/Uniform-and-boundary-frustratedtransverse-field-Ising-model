# Clean XY model with PBC by exact diagonalization method
# Creator: Viet-Anh Tran - anhvt.tran@gmail.com

import numpy as np
import math as mt 
import matplotlib.pyplot as plt
from numpy import linalg as LA
import time
import numba 
from datetime import datetime
import winsound
import json

# Initialization
PI = mt.pi
sigma_x = np.matrix([[0,1],[1,0]]) #Pauli matrices
sigma_y = np.matrix([[0,-1j],[1j,0]])
sigma_z = np.matrix([[1,0],[0,-1]])
I = np.matrix([[1,0],[0,1]]) #Identity matrix

def write_dict(a_dict,dict_name): # to save into json file
    print("Start writing dictionary data into a json file")
    with open(dict_name, "w") as fp:
        json.dump(a_dict,fp)
        print("Done writing JSON data into .json file")

# Hamiltonian
#@numba.njit
def Hamiltonian(L,h,epsilon):
    # Transverse Part.
    i=1
    while i<=L:
        if i==1:
            t=sigma_z
        else :
            t = I
        j=2
        while j<=L:
            if i==j:
                t=np.kron(sigma_z,t)
            else:
                t=np.kron(I,t)
            j+=1
        if i==1:
            transverse = t
        else:
            transverse = transverse + t
        i+=1
        
  # Interaction x- part.    
    i=1
    while i<=L:
        if i==1 or i==L:
           t=sigma_x
        else:
           t = I
        j=2
        while j<=L:
            if i==j or j==i+1:
                t=np.kron(sigma_x,t)
            else:
                t=np.kron(I,t)
            j+=1
        if i==1:
            interaction_x = t
        else:
            interaction_x= interaction_x + t
        i+=1
        
    # Interaction y- part.    
    i=1
    while i<=L:
        if i==1 or i==L:
           t = sigma_y
        else:
           t = I
        j=2
        while j<=L:
            if i==j or j==i+1:
                t = np.kron(sigma_y,t)
            else:
                t = np.kron(I,t)
            j+=1
        if i==1:
            interaction_y = t
        else:
            interaction_y = interaction_y + t
        i+=1

    # Additional energy term 
    i=1
    while i<=L:
        if i==1:
            t=sigma_x
        else :
            t = I
        j=2
        while j<=L:
            if i==j:
                t=np.kron(sigma_x,t)
            else:
                t=np.kron(I,t)
            j+=1
        if i==1:
            Additional_Energy = t
        else:
            Additional_Energy = Additional_Energy + t
        i+=1

    #Total Hamiltonian
    #epsilon = -10**-3
    H_xy = -( h*transverse + interaction_x) + epsilon*Additional_Energy# + interaction_y )  
    return(H_xy)

# Diagonalization function
def energy_levels(A):
    EigenEnergies,EigenVectors = LA.eigh(A)
    #EigenEnergiesSorted=np.sort(EigenEnergies)
    return [EigenEnergies,EigenVectors]

# Function of magnetization matrix 
def Magnetization_matrix(L):
    i=1
    while i<=L:
        if i==1:
            t=sigma_x
        else :
            t = I
        j=2
        while j<=L:
            if i==j:
                t=np.kron(sigma_x,t)
            else:
                t=np.kron(I,t)
            j+=1
        if i==1:
            Magnetization = t
        else:
            Magnetization = Magnetization + t
        i+=1
    return Magnetization/L


################################################################

# Execution
if __name__== '__main__':
    epsilon_values=[-10**-3,0,10**-3]
    h_values=np.arange(0,2.05,0.05).tolist() # g runs from 0 to 2.05 with the step size is 0.05
    L_values=[4,6,8,10,12] # Values of L (Number of sites) => dimension of the Hamiltonian = 2**L
   
    ## preallocate lists
    Egs_values=[]
    First_excited_values=[]
    Second_excited_values=[]
    Mag_values=[]
    
    start_time = datetime.now() # Start counting the time 
    for epsilon in epsilon_values:
        Egs_epsilon_const=[]
        First_excited_epsilon_const=[]
        Second_excited_epsilon_const=[]
        Mag_epsilon_const=[]

        for L in L_values:
            Egs_L_const=[]
            First_excited_L_const=[]
            Second_excited_L_const=[]
            Mag_L_const=[]
            for h in h_values:
                Es,Ev = energy_levels( Hamiltonian(L,h,epsilon/L) )
                Egs_L_const.append( Es[0] ) # GS energy at one N value   
                First_excited_L_const.append( Es[1] )
                Second_excited_L_const.append( Es[2] ) 
                Mag_L_const.append( np.around (np.dot( np.dot( Ev[:,0].getH(), Magnetization_matrix(L) ) , Ev[:,0] ).item(),6) )
            Egs_epsilon_const.append( Egs_L_const ) 
            First_excited_epsilon_const.append( First_excited_L_const )
            Second_excited_epsilon_const.append( Second_excited_L_const )
            Mag_epsilon_const.append( Mag_L_const )
        Egs_values.append( Egs_epsilon_const ) 
        First_excited_values.append( First_excited_epsilon_const )
        Second_excited_values.append( Second_excited_epsilon_const )
        Mag_values.append( Mag_epsilon_const )
    
    # Writing out the data 
    d = {"epsilon_values" : epsilon_values,"epsilon_values" : epsilon_values,"L_values" : L_values , "h_values" : h_values, "Egs_values" : Egs_values, "First_excited_values": First_excited_values , "Second_excited_values":Second_excited_values,"Mag_values":Mag_values }
    write_dict(d,"ed.json")
    end_time = datetime.now()
    print("Calculation completed")
    print('Duration: {}'.format(end_time - start_time))
    
winsound.Beep(500, 2500) # alarm when the code finished already

# or 
# import os 
# os.system('play -nq -t alsa synth {} sine {}'.format(5, 500))
# for Linux and Mac 