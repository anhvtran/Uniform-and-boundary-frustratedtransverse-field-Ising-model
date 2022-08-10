# Definition here is a bit different from the conventional Ising model for the convenient purpose when running J_link
# positive J for ANTIFERROMAGNETICS case and negative J for FERROMAGNETICS case (in case of ferro, you have to change the sign of J_link to J as well)
# keep positive J in this study currently in order to observe the antiferro

from __future__ import absolute_import, division, print_function
import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.fftpack import hilbert
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import time
import numba 
from datetime import datetime
import winsound
import json

sigma_x = np.matrix([[0,1],[1,0]])
sigma_z = np.matrix([[1,0],[0,-1]])
I = np.matrix([[1,0],[0,1]]) #Identity matrix

def write_dict(a_dict,dict_name): # to save into json file
    print("Start writing dictionary data into a json file")
    with open(dict_name, "w") as fp:
        json.dump(a_dict,fp)
        print("Done writing JSON data into .json file")

def hilbertspace_dimension(L):
        ''' return dimension of hilbertspace '''
        return 2**L

def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> (L-1-site))  & 1

def get_hamiltonian_sparse(L, J, J_link, hz, epsilon):
    '''
    Creates the Hamiltonian of the Transverse Field Ising model
    '''
    # Initialization
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    for state in range(hilbertspace_dimension(L)):                    
        for site in range(L-1):
            if ( (get_site_value( state, site) == 0) and (get_site_value( state, site + 1) == 0) ):
                new_state = state + 2**(L-site-1) + 2**(L-site-2)
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append(J)
            if ( (get_site_value( state, site) == 0) and (get_site_value( state, site + 1) == 1) ): 
                new_state = state + 2**(L-site-1) - 2**(L-site-2)
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append(J) 
            if ( (get_site_value( state, site) == 1) and (get_site_value( state, site + 1) == 1) ):
                new_state = state - 2**(L-site-1) - 2**(L-site-2)
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append(J)
            if ( (get_site_value( state, site) == 1) and (get_site_value( state, site + 1) == 0) ):
                new_state = state - 2**(L-site-1) + 2**(L-site-2)
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append(J)

        # The last bond 
        if ( (get_site_value( state, L-1) == 0) and (get_site_value( state, 0) == 0) ):
            new_state = state + 2**0 + 2**(L-1)
            hamiltonian_rows.append(new_state) 
            hamiltonian_cols.append(state)
            hamiltonian_data.append(J_link)
        if ( (get_site_value( state, L-1) == 0) and (get_site_value( state, 0) == 1) ): 
            new_state = state + 2**0 - 2**(L-1)
            hamiltonian_rows.append(new_state) 
            hamiltonian_cols.append(state)
            hamiltonian_data.append(J_link)
        if ( (get_site_value( state, L-1) == 1) and (get_site_value( state, 0) == 1) ):
            new_state = state - 2**0 - 2**(L-1)
            hamiltonian_rows.append(new_state)
            hamiltonian_cols.append(state)
            hamiltonian_data.append(J_link)
        if ( (get_site_value( state, L-1) == 1) and (get_site_value( state, 0) == 0) ):
            new_state = state - 2**0 + 2**(L-1)
            hamiltonian_rows.append(new_state)
            hamiltonian_cols.append(state)
            hamiltonian_data.append(J_link)

        # # Apply the perturbation term used for negative J case (ferro)
        # for site in range(L):
        #     if get_site_value( state,site ) == 1:
        #         new_state =  state - 2**(L-site-1)
        #         hamiltonian_rows.append(new_state)
        #         hamiltonian_cols.append(state)
        #         hamiltonian_data.append(epsilon)
        #     else:
        #         new_state =  state + 2**(L-site-1)
        #         hamiltonian_rows.append(new_state)
        #         hamiltonian_cols.append(state)
        #         hamiltonian_data.append(epsilon) 

        # Apply the perturbation term used for positive J case (antiferro)
        for site in range(L):
            if get_site_value( state,site ) == 1:
                new_state =  state - 2**(L-site-1)
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append( (-1)**(site + 1) * epsilon )
            else:
                new_state =  state + 2**(L-site-1)
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append( (-1)**(site + 1) * epsilon ) 

        # Apply transverse field
        diagonal_value = 0
        for site in range(L):
            if get_site_value(state, site) == 1:
                diagonal_value += -hz
            else:
                diagonal_value -= -hz

        hamiltonian_rows.append(state)
        hamiltonian_cols.append(state)
        hamiltonian_data.append(diagonal_value)

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data

# Function of magnetization matrix 
def Magnetization_matrix(L): # for negative J case
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
            Magnetization += t
        i+=1
    return Magnetization/L

def Staggered_Magnetization_matrix(L): # for positive J case
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
            Staggered_Magnetization = -t
        else:
            Staggered_Magnetization += (-1)**(i) * t 
        i+=1
    return Staggered_Magnetization/L     

def magnetization_sparse(L): # when J negative 
    magnetization_rows = []
    magnetization_cols = []
    magnetization_data = []
    for state in range( hilbertspace_dimension(L) ):
        for site in range( L ):
            if get_site_value(state,site) == 0:
                new_state = state + 2**(L-site-1)
                magnetization_rows.append(new_state)
                magnetization_cols.append(state)
                magnetization_data.append(1/L)
            else:
                new_state = state - 2**(L-site-1)
                magnetization_rows.append(new_state)
                magnetization_cols.append(state)
                magnetization_data.append(1/L)
    return magnetization_rows,magnetization_cols,magnetization_data

def staggered_magnetization_sparse(L): # when J postitive
    smagnetization_rows = []
    smagnetization_cols = []
    smagnetization_data = []
    for state in range( hilbertspace_dimension(L) ):
        for site in range( L ):
            if get_site_value(state,site) == 0:
                new_state = state + 2**(L-site-1)
                smagnetization_rows.append(new_state)
                smagnetization_cols.append(state)
                smagnetization_data.append( (-1)**(site+1) * 1/L )
            else:
                new_state = state - 2**(L-site-1)
                smagnetization_rows.append(new_state)
                smagnetization_cols.append(state)
                smagnetization_data.append( (-1)**(site+1) * 1/L )
    return smagnetization_rows,smagnetization_cols,smagnetization_data

def get_correlation_sparse( L,site1,site2 ):
    '''
    Create the sparse form of correlation function between two spins
    '''
    correlation_rows = []
    correlation_cols = []
    correlation_data = []

    for state in range( hilbertspace_dimension(L) ):
        if ( (get_site_value( state, site1) == 0) and (get_site_value( state, site2) == 0) ):
            new_state = state + 2**(L-site1-1) + 2**(L-site2-1)
            correlation_rows.append(new_state)
            correlation_cols.append(state)
            correlation_data.append(1)
        if ( (get_site_value( state, site1) == 0) and (get_site_value( state, site2) == 1) ): 
            new_state = state + 2**(L-site1-1) - 2**(L-site2-1)
            correlation_rows.append(new_state)
            correlation_cols.append(state)
            correlation_data.append(1)
        if ( (get_site_value( state, site1) == 1) and (get_site_value( state, site2) == 1) ):
            new_state = state - 2**(L-site1-1) - 2**(L-site2-1)
            correlation_rows.append(new_state)
            correlation_cols.append(state)
            correlation_data.append(1)
        if ( (get_site_value( state, site1) == 1) and (get_site_value( state, site2) == 0) ):
            new_state = state - 2**(L-site1-1) + 2**(L-site2-1)
            correlation_rows.append(new_state)
            correlation_cols.append(state)
            correlation_data.append(1)
    return correlation_rows, correlation_cols, correlation_data
    
# Execution
if __name__== '__main__':

    epsilon_values=[-10**-3, 0, 10**-3]
    L_values = [8,10,12,14,16]
    J = 1 # antiferro 
    #J = -1 # ferro
    n_lowest_eigenvalues = 3
    hz_values = [0.2 , 0.5] # just in case
    J_link_values = np.linspace(0,3,150)

    start_time = datetime.now() # Start counting the time
    
    hz = 0.2
    #hz = 0.7 
    es0_values = []
    es1_values = []
    es2_values = []
    mag_values = [] 
    corr_values = []
    for epsilon in epsilon_values:
        es0_epsilon_values = [] 
        es1_epsilon_values = []
        es2_epsilon_values = []
        mag_epsilon_values = [] 
        corr_epsilon_values = []
        for L in L_values: # run L 
            es0_L_const = []
            es1_L_const = []
            es2_L_const = []
            mag_L_const = []
            corr_L_const = []
            for J_link in J_link_values: # run J_link
                rows, cols, data = get_hamiltonian_sparse(L, J, J_link, hz, epsilon/L)
                hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)),shape=(2**L, 2**L))
                es,ev = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                                            which='SA', return_eigenvectors=True,
                                            maxiter=1000000) 
                es0_L_const.append( es[0] )
                es1_L_const.append( es[1] )
                es2_L_const.append( es[2] )
                #mag_L_const.append( np.around (np.dot( np.dot( np.ndarray.conjugate(ev[:,0]), Magnetization_matrix(L) ) , ev[:,0] ).item(),6) ) # negative J case 
                #mag_L_const.append( np.around (np.dot( np.dot( np.ndarray.conjugate(ev[:,0]), Staggered_Magnetization_matrix(L) ) , ev[:,0] ).item(),6) ) # positive J case
                # magnetization sparse 
                # mag_rows, mag_cols, mag_data = magnetization_sparse(L)
                # mag_matrix = sp.sparse.csr_matrix((mag_data, (mag_rows, mag_cols)),shape=(2**L, 2**L))
                # mag_L_const.append( np.around( np.dot(ev[:,0], mag_matrix.dot(ev[:,0])).item(),6) )
                # staggered magnetization sparse
                smag_rows, smag_cols, smag_data = staggered_magnetization_sparse(L) 
                smag_matrix = sp.sparse.csr_matrix((smag_data, (smag_rows, smag_cols)),shape=(2**L, 2**L))
                mag_L_const.append( np.around( np.dot(ev[:,0], smag_matrix.dot(ev[:,0])).item(),6) )
                # correlation sparse
                corr_rows1, corr_cols1, corr_data1 = get_correlation_sparse(L,0,L-1)  
                corr_matrix1 = sp.sparse.csr_matrix( (corr_data1, (corr_rows1, corr_cols1)),shape=(2**L, 2**L) )
                corr1 = np.dot(ev[:,0], corr_matrix1.dot(ev[:,0])) 
                #print( 'corr1' ) ; print(corr1)
                corr_rows2, corr_cols2, corr_data2 = get_correlation_sparse(L,L//2, (L//2)+1 ) 
                corr_matrix2 = sp.sparse.csr_matrix( (corr_data2, (corr_rows2, corr_cols2)),shape=(2**L, 2**L) )
                corr2 = np.dot(ev[:,0], corr_matrix2.dot(ev[:,0]))
                #print( 'corr2' ) ; print(corr2)
                corr_L_const.append( corr1/corr2 )

            es0_epsilon_values.append( es0_L_const )
            es1_epsilon_values.append( es1_L_const )
            es2_epsilon_values.append( es2_L_const )
            mag_epsilon_values.append( mag_L_const )
            corr_epsilon_values.append( corr_L_const )

        es0_values.append( es0_epsilon_values )
        es1_values.append( es1_epsilon_values )
        es2_values.append( es2_epsilon_values )
        mag_values.append( mag_epsilon_values )
        corr_values.append( corr_epsilon_values )

    # Writing out the data 
    d = {"epsilon_values" : epsilon_values,"epsilon_values" : epsilon_values,"L_values" : L_values , "hz_values" : hz_values, "J_link_values" : J_link_values.tolist(), "es0_values" : es0_values, "es1_values": es1_values , "es2_values":es2_values,"mag_values":mag_values, "corr_values":corr_values }
    write_dict(d,"full_evenL_upto16_h0.2.json")

    # Plotting
    color = 'cgbrmykcgbrmyk' # color codes 

    # 1st gap at a given hz
    es0_h = es0_values[1]
    es1_h = es1_values[1]
    es2_h = es2_values[1]
    for i in  range(0,len(L_values)):
        C = str(color[i])
        L_lab = 'L = '+ str(L_values[i])
        plt.plot(J_link_values,np.array(es1_h[i])-np.array(es0_h[i]), color=C, label = str(L_lab) )
        plt.legend(loc= "upper right", prop={'size': 6})
        #plt.title("Energy gap between First Excited State and Ground State when h = %s" %hz)
    plt.xlabel("$J_{link}$")
    plt.ylabel("$E_1 - E_0$")
    #plt.savefig('10.png',dpi=300)
    plt.show()

    # 2nd gap at a given hz 
    for i in  range(0,len(L_values)):
        C = str(color[i])
        L_lab = 'L = '+ str(L_values[i])
        plt.plot(J_link_values,np.array(es2_h[i])-np.array(es0_h[i]), color=C, label = str(L_lab) )
        plt.legend(loc= "upper right", prop={'size': 6})
        #plt.title("Energy gap between Second Excited State and Ground State when h = %s" %hz)
    plt.xlabel("$J_{link}$")
    plt.ylabel("$E_2 - E_0$")
    #plt.savefig('11.png',dpi=300)
    plt.show()
    
    # Plot the magnetization or staggered magnetization at a given hz 
    line_values = ['solid','dotted','dashed']
    for j in range(0,len(epsilon_values)):
        line = str(line_values[j])
        M = mag_values[j]
        for i in  range(0,len(L_values)):
            C = str(color[i])
            L_lab = 'L = '+ str(L_values[i])
            plt.plot(J_link_values,M[i], color=C,ls = line, label = str(L_lab) if j == 0 else "" )
            plt.legend(loc= "upper right", prop={'size': 6})
            #plt.title("GS Magnetization when h = %s" %hz) 
        plt.xlabel("$J_{link}$")
        #plt.ylabel("$M$")
        plt.ylabel("$M_{stag}$")
    #plt.savefig('12.png',dpi=300)
    plt.show()

    # Plot the correlation relation 
    line_values = ['solid','dotted','dashed']
    for j in range(0,len(epsilon_values)):
        line = str(line_values[j])
        corr = corr_values[j]
        for i in  range(0,len(L_values)):
            C = str(color[i])
            L_lab = 'L = '+ str(L_values[i])
            plt.plot(J_link_values,corr[i], color=C,ls = line, label = str(L_lab) if j == 0 else "" )
            plt.legend(loc= "upper right", prop={'size': 6})
            #plt.title("GS Magnetization when h = %s" %hz) 
        plt.xlabel("$J_{link}$")
        #plt.ylabel("$M$")
        plt.ylabel("$C$")
    plt.savefig('corr_13.png',dpi=300)
    plt.show()

    end_time = datetime.now()
    print("Calculation completed")
    print('Duration: {}'.format(end_time - start_time))
    print("Run successful!")
winsound.Beep(500, 500) # alarm when the code finished already

# or 
# import os 
# os.system('play -nq -t alsa synth {} sine {}'.format(5, 500))
# for Linux and Mac 