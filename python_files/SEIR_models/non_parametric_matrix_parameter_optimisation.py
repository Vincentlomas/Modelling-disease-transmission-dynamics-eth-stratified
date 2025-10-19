# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 10:10:16 2024

@author: Vincent Lomas

Optimisation of model parameters
"""
import scipy
import numpy as np
from thesis_modules import *


def model_simulation_non_parametric(a0, CAR, SEIR_0, N_vec, time, sigma,gamma, eth_data=None):
    
    # Convert a to a column vector
    a0.flatten()
    a0=np.array([a0]).T
    
    # beta is the transmission matrix
    beta = non_parametric_per_capita_transmission_matrix(a0, N_vec, eth_data)
     
    solution = scipy.integrate.solve_ivp(SEIR_model,[0,time],SEIR_0,args=(beta,sigma,gamma))
    
    attack_rates = return_attack_rates(CAR)
    attack_rate_percent_diff = attack_rates.flatten() - 100 * (N_vec.flatten() - solution.y[0:4,-1] - solution.y[4:8,-1] - solution.y[8:12,-1]) / N_vec.flatten()
    
    return attack_rate_percent_diff



def non_parametric_optimisation(N_vec, N_vec_vacc, N_vec_vacc_boosted, is_vacc, time, sigma, gamma):

    ### initial group population distributions
    # 0.01% initial expoure for each group
    S, Sv, Svb, E, Ev, Evb, I, Iv, Ivb, R, In = initial_group_populations(N_vec,is_vacc,N_vec_vacc, N_vec_vacc_boosted)
    
    # loop over both statisticals areas (if required)
    # for is_SA1 in [False]: #(False,True):
        
    #     if is_SA1:
    #         is_statsnz_list = [True] # [True]
    #         print_statment0 = 'SA1 '
    #     else:
    #         is_statsnz_list = [True,False] # [True, False]
    print_statment0 = 'SA2 '
    is_SA1 = False
    is_statsnz = False
    
    ### import data
    eth_data = SA_import(is_SA1, is_statsnz)
    
    # Append print statement
    if is_statsnz:
        print_statment1 = print_statment0 + 'total '
    else:
        print_statment1 = print_statment0 + 'priority '
    
    # Loop over each case asceratinment rate
    for CAR in ['4060',40,50,60]:
        # Append case ascertainment rate to print statement
        print_statement = print_statment1 + str(CAR)
      
        if True: # Skip specific cases if wanted here
            # initial transmission vector
            a = np.array([[1,1,1,1]]).T
    
                
            #     # Runs if transmission rates are negaive or error is too large
            sol = scipy.optimize.least_squares(model_simulation_non_parametric,
                                               (0.3,0.5,0.4,1),
                                               bounds = (0.1,2), args = (CAR, np.concatenate([S,Sv,Svb,E, Ev, Evb, I, Iv, Ivb,R,In]), N_vec, time,sigma,gamma, eth_data))
            result = sol.x
                
            
            a = np.array([result]).T
            
            
            err = np.sum(abs(np.array(model_simulation_non_parametric(a.flatten(), CAR, np.concatenate([S,Sv,Svb,E, Ev, Evb, I, Iv, Ivb,R,In]),
                                                       N_vec, time, sigma,gamma, eth_data))))
            
            
            print(print_statement)
            print()
            print(f'Transmission rates {a}')
            print(f' Error: {err}')
            
            # print results
            print()
            if err > 10**-5:
                print('-'*97)
                print(f'Warning: Above set up has error of {err}')
                print('-'*97)
            
            save_epsilon_a_vector(0,a,CAR,is_vacc,is_prop = False,is_non_parametric = True,
                                      is_SA1=is_SA1,is_statsnz=is_statsnz, is_old = False)
