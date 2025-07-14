# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 10:10:16 2024

@author: Vincent Lomas

Optimisation of model parameters
"""
import scipy
import numpy as np
from thesis_modules import SEIR_model, save_epsilon_a_vector, initial_group_populations, return_attack_rates

def model_simulation_prop(a0, CAR, SEIR_0, N_vec, time, sigma,gamma):
    '''Acts as a loss function for fitting a contact rate to a value of
    assortativity.
    Takes as an input a contact vector (of floats). Also uses a global float 
    variable, epsilon, as an assortativity value.
    It then runs an SEIR model and determines the difference between some
    specified global variable attack rates (mao_N2,pac_N2,asi_N2,oth_N2)
    the absolute values of these differences is added and returned'''
    # Conver a to a column vector
    a0.flatten()
    a0=np.array([a0]).T
    
    # beta is the transmission matrix
    beta = (a0 @ a0.T)/(a0.T @ N_vec)
     
    solution = scipy.integrate.solve_ivp(SEIR_model,[0,time],SEIR_0,args=(beta,sigma,gamma))
    
    attack_rates = return_attack_rates(CAR)
    attack_rate_percent_diff = attack_rates.flatten() - 100 * (N_vec.flatten() - solution.y[0:4,-1] - solution.y[4:8,-1]) / N_vec.flatten()
    
    return attack_rate_percent_diff


def proportionate_optimisation(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma):
    ### initial group population distributions
    # 0.01% initial expoure for each group
    S, Sv, E, I, R, In = initial_group_populations(N_vec,is_vacc,N_vec_vacc)
    
    # Loop over each case asceratinment rate
    for CAR in ['4060', 40, 50, 60]:
        
        # initial guess for transmission vector - convergence can sometime be sensitive
        # a = a_vector_prop_mixing(CAR,is_vacc)
        a=np.ones([1,4])
    
        if True: # Skip specific cases if wanted here
        
            rates_all_positive = False
            
            while not rates_all_positive:
                result = scipy.optimize.fsolve(model_simulation_prop,a.flatten(), args=(CAR, np.concatenate([S,Sv,E,I,R,In]), N_vec, time, sigma,gamma))
                
                a = np.array([result])
                
                if np.all(a>0):
                    rates_all_positive = True
                else:
                    a=np.random.uniform(0,1,[4,1])
            
            print(CAR)
            print()
            print(f'Transmission rates {a}')
            print()
            
            # save results
            save_epsilon_a_vector(0,a,CAR,is_vacc,is_prop = True,is_non_parametric = False)