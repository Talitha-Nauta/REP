import numpy as np
import time
import os
import gurobipy as gp
from scipy.stats import binom,geom,chi2

#global parameters
significance_level = 0.01 # number of false alarms, type I error, same for all tests
a_min = 0 #No minimum attack length

#define closed loop systems
def define_system(system):    
    if system == 'furuta_zero':
        # | Furuta Pendulumn 3 states zero
        Phi_miss = np.array([[1.00259871763614,	0.00500433044615377, 0,	0,	0,	0,	-0.0843405496542633],[1.03993704285104,	1.00259871763614,	0,	0,	0,	0,	-33.7508240911965],[-0.0675434992042660,	-0.000168785681590098,	1,	0,	0,	0,	39.2130808215931],[0,	0,	0,	1,	0,	0,	0],[0,	0,	0,	0,	1,	0,	0],[0,	0,	0,	0,	0,	1,	0],[0,	0,	0,	0,	0,	0,	0]])
        Phi_hit = np.array([[1.00259871763614,	0.00500433044615377,	0,	0,	0,	0,	-0.0843405496542633],[1.03993704285104,	1.00259871763614,	0,	0,	0,	0,	-33.7508240911965],[-0.0675434992042660,	-0.000168785681590098,	1,	0,	0,	0,	39.2130808215931],[1.00259871763614,	0.00500433044615377,	0,	0,	0,	0,	-0.0843405496542633],[1.03993704285104,	1.00259871763614,	0,	0,	0,	0,	-33.7508240911965],[-0.0675434992042660,	-0.000168785681590098,	1,	0,	0,	0,	39.2130808215931],[0,	0,	0,	0.428000000000000,	0.0307000000000000,	0.0119000000000000,	0]])
        C = np.eye(3) 
        nx = Phi_hit.shape[1]  # number of states (closed-loop system)
        Q = np.zeros((nx, nx)) 
        Q[0, 0] = 1
        x0 = np.array([0.0097627,  0.04303787, 0.02055268, 0.0097627,  0.04303787, 0.02055268,  0])
    elif system == 'furuta_hold':
        # | Furuta Pendulumn 3 states hold
        Phi_miss = np.array([[1.00259871763614,	0.00500433044615377, 0,	0,	0,	0,	-0.0843405496542633],[1.03993704285104,	1.00259871763614,	0,	0,	0,	0,	-33.7508240911965],[-0.0675434992042660,	-0.000168785681590098,	1,	0,	0,	0,	39.2130808215931],[0,	0,	0,	1,	0,	0,	0],[0,	0,	0,	0,	1,	0,	0],[0,	0,	0,	0,	0,	1,	0],[0,	0,	0,	0,	0,	0,	1]])
        Phi_hit = np.array([[1.00259871763614,	0.00500433044615377,	0,	0,	0,	0,	-0.0843405496542633],[1.03993704285104,	1.00259871763614,	0,	0,	0,	0,	-33.7508240911965],[-0.0675434992042660,	-0.000168785681590098,	1,	0,	0,	0,	39.2130808215931],[1.00259871763614,	0.00500433044615377,	0,	0,	0,	0,	-0.0843405496542633],[1.03993704285104,	1.00259871763614,	0,	0,	0,	0,	-33.7508240911965],[-0.0675434992042660,	-0.000168785681590098,	1,	0,	0,	0,	39.2130808215931],[0,	0,	0,	0.428000000000000,	0.0307000000000000,	0.0119000000000000,	0]])
        C = np.eye(3)
        nx = Phi_hit.shape[1]  # number of states (closed-loop system)
        Q = np.zeros((nx, nx)) 
        Q[0, 0] = 1
        x0 = np.array([0.0097627,  0.04303787, 0.02055268, 0.0097627,  0.04303787, 0.02055268,  0])
    elif system == 'qt_nmp':
        # | quad tank process non-minimum phase LTH!, with Ts=0.1 and hold:
        Phi_miss = np.array([[0.995721279601902,		0,		0.00609376139131610,		0,		0,		0,		0,		0,		0,		0,		0.00977494639256366,		6.97640108027585e-05],[0,		0.995721279601902,		0,		0.00609376139131610,		0,		0,		0,		0,		0,		0,		6.97640108027585e-05,		0.00977494639256366],[0,		0,		0.993893151184507,		0,		0,		0,		0,		0,		0,		0,		0,		0.0227872790460744],[0,		0,		0,		0.993893151184507,		0,		0,		0,		0,		0,		0,		0.0227872790460744,		0],[0,		0,		0,		0,		1,		0,		0,		0,		0,		0,		0,		0],[0,		0,		0,		0,		0,		1,		0,		0,		0,		0,		0,		0],[0,		0,		0,		0,		0,		0,		1,		0,		0,		0,		0,		0],[0,		0,		0,		0,		0,		0,		0,		1,		0,		0,		0,		0],[0,		0,		0,		0,		0,		0,		0,		0,		1,		0,		0,		0],[0,		0,		0,		0,		0,		0,		0,		0,		0,		1,		0,		0],[0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		1,		0],[0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		1]])
        Phi_hit = np.array([[0.995721279601902,	0,	0.00609376139131610,	0,	0,	0,	0,	0,	0,	0,	0.00977494639256366,	6.97640108027585e-05],[0,	0.995721279601902,	0,	0.00609376139131610,	0,	0,	0,	0,	0,	0,	6.97640108027585e-05,	0.00977494639256366],[0,	0,	0.993893151184507,	0,	0,	0,	0,	0,	0,	0,	0,	0.0227872790460744],[0,	0,	0,	0.993893151184507,	0,	0,	0,	0,	0,	0,	0.0227872790460744,	0],[0.497860639800951,	0,	0.00304688069565805,	0,	0,	0,	0,	0,	0,	0,	0.00488747319628183,	3.48820054013793e-05],[0,	0.497860639800951,	0,	0.00304688069565805,	0,	0,	0,	0,	0,	0,	3.48820054013793e-05,	0.00488747319628183],[0,	0,	0,	0,	0.772881693607458,	-0.0169763980536733,	0.599775679106704,	-0.00925928466940067,	0.0117654040072148,	-0.00982881503410122,	0,	0],[0,	0,	0,	0,	-0.0169763980536730,	0.772881693607459,	-0.00925928466940027,	0.599775679106705,	-0.00982881503410170,	0.0117654040072148,	0,	0],[0,	0,	0,	0,	0.277581805155136,	-0.0127803678293101,	-0.180007660158698,	-0.0154730590871807,	0.970884758508424,	0.0133859005910090,	0,	0],[0,	0,	0,	0,	-0.0127803678293110,	0.277581805155135,	-0.0154730590871815,	-0.180007660158699,	0.0133859005910091,	0.970884758508425,	0,	0],[0,	0,	0,	0,	-0.560855370378794,	-1.73272260054955,	-0.679021793514529,	-0.942400399522460,	0.587428651044457,	-1.00970338009905,	0,	0],[0,	0,	0,	0,	-1.73272260054953,	-0.560855370378755,	-0.942400399522418,	-0.679021793514497,	-1.00970338009910,	0.587428651044452,	0,	0]]) 
        C = np.array([[0.5, 0, 0, 0],[0, 0.5, 0, 0]])       
        nx = Phi_hit.shape[1]
        Q = np.zeros((nx, nx)) 
        Q[0, 0] = 1
        Q[1, 1] = 1
        x0 = np.array([0.09762701, 0.43037873, 0.20552675, 0.08976637, 0.048813505,  0.215189365, 0, 0, 0, 0, 0, 0])
    else:
        print("This system does not exist")
        Phi_miss = -1, 
        Phi_hit = -1
        C = -1
        Q = -1
        x0 = -1
        nx = -1
        
    return Phi_miss, Phi_hit, C, Q, x0, nx

#Create output folders
try:
    os.makedirs("results_bin_chi2/furuta_zero")
    os.makedirs("results_bin_chi2/qt_nmp")
    os.makedirs("results_bin/furuta_zero")
    os.makedirs("results_bin/qt_nmp")
    os.makedirs("results_chi2/furuta_zero")
    os.makedirs("results_chi2/qt_nmp")
    print("starting")

except OSError:  # parts of this experiment set might already have been performed
    print("part of the results might already present - remove existing result folders")
    exit(1)

window_lengths = [20] #Last step not included in arange
prob = 0.01*np.arange(9,26,8) #Avoid numerical error #steps of two is enough to cover all combinations for allowed deadline misses per window_length

for a_length in np.arange(1,28,1): #Corresponds to lengths 1-25 in pratice as both first and last are forced to be a deadline hit in the optimisation implementation
    st = time.time() #Start time
    for system in ['furuta_zero','qt_nmp']:  #furuta_zero, qt_nmp
        for a_max in [4]:  #a_max 4
            for window_length in window_lengths: #20
                for p_miss in prob: #p_miss 0.09:0.08:0.25

                    Phi_miss, Phi_hit, C, Q, x0, nx = define_system(system)
                    #Do bin + chi2 test
                    try:
                        exec(open("./optimise_bin_chi2.py").read()) 
                    except:
                        print(f"There was a problem for length={a_length}, system={system}, a_max = {a_max}, window length = {window_length}, p miss = {p_miss}, test = bin + chi2")

                    #Do binary test
                    try:
                        exec(open("./optimise_bin.py").read()) 
                    except:
                        print(f"There was a problem for length={a_length}, system={system}, a_max = {a_max}, window length = {window_length}, p miss = {p_miss}, test = bin")

                    #Do chi2 test
                    try:
                        exec(open("./optimise_chi2.py").read()) 
                    except:
                        print(f"There was a problem for length={a_length}, system={system}, a_max = {a_max}, window length = {window_length}, p miss = {p_miss}, test = chi2")

    #Print elapsed time per window length
    et = time.time() 
    elapsed_time = et - st
    print('Execution time for sequence length', a_length, 'is: ', elapsed_time, 'seconds')
