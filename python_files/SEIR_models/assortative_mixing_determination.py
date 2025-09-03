# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 20:44:38 2024

@author: Vincent Lomas

Fitting assortative transmission matrix using census data (SA1 2018 data)

Data from: https://www.stats.govt.nz/information-releases/statistical-area-1-dataset-for-2018-census-updated-march-2020
Part 1 data file was used - Ethnic group (grouped total responses)
"""

import os
import scipy
import numpy as np
import matplotlib.pyplot as plt
from thesis_modules import *

import warnings # Supress warnings of duplicate function calls during minimisation

warnings.filterwarnings("ignore") # Ignore all warnings --- BAD IDEA

### SOme error

def epsilon_variance(e0,a, N_vector_eth, C_prime):
    C_temp = (1-e0)*((a @ a.T)*(N_vector_eth@N_vector_eth.T))/(a.T @ N_vector_eth) + e0 * np.diag(a[:,0] * N_vector_eth[:,0])
    return np.sum(np.abs(C_temp-C_prime))


def model_simulation_assortative(a0, CAR, SEIR_0, N_vec, time, sigma,gamma, epsilon):
    '''Acts as a loss function for fitting a contact rate to a value of
    assortativity.
    Takes as an input a contact vector (of floats). Also uses a global float 
    variable, epsilon, as an assortativity value.
    It then runs an SEIR model and determines the difference between some
    specified global variable attack rates (mao_N,pac_N,asi_N,oth_N)
    the absolute values of these differences is added and returned'''
    # Conver a to a column vector
    a0.flatten()
    a0=np.array([a0]).T 
    
    # beta is the transmission matrix
    beta = (1-epsilon)*(a0 @ a0.T)/(a0.T @ N_vec) + epsilon * np.diag(a0[:,0] / N_vec[:,0])
    
    solution = scipy.integrate.solve_ivp(SEIR_model,[0,time],SEIR_0,args=(beta,sigma,gamma))
    
    
    attack_rates = return_attack_rates(CAR)
    attack_rate_percent_diff = attack_rates.flatten() - 100 * (N_vec.flatten() - solution.y[0:4,-1] - solution.y[4:8,-1]) / N_vec.flatten()
    
    return attack_rate_percent_diff


### initialization of parameters being determined
#a = np.array([[0.5,0.5,0.5,0.5]]).T

### Old populations
# mao_N =  802030 # From dataset cases/cases_per_100k
# pac_N =  390970
# asi_N =  834055
# oth_N = 3167397
# tot_N = mao_N + asi_N + pac_N + oth_N
# N_vec = np.array([[mao_N],
#          [pac_N],
#          [asi_N],
#          [oth_N]])

def assortative_optimisation(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma):
    tol = 10**-6
    max_count = 30
    
    for is_SA1, is_statsnz in [(True,True),(False,True),(False,False)]:
            
            if is_SA1:
                print_statement0 = "SA1 "
            else:
                print_statement0 = "SA2 "
            
            if is_statsnz:
                print_statement1 = print_statement0 + "total"
            else:
                print_statement1 = print_statement0 + "priority"
                
            # Import data
            eth_data = SA_import(is_SA1, is_statsnz)
        
            # Draw populations from vectors
            mao_N_eth = np.sum(eth_data[:,0])
            pac_N_eth = np.sum(eth_data[:,1])
            asi_N_eth = np.sum(eth_data[:,2])
            oth_N_eth = np.sum(eth_data[:,3])
            N_vector_eth = np.array([[mao_N_eth,pac_N_eth,
                                  asi_N_eth,oth_N_eth]],dtype='int64').T
            
            ### iterate over all case ascertainment rates
            for CAR in ['4060',40,50,60]:
                print_statement = print_statement1 + f' CAR{CAR}'
                
                if True: # Skip cases if wanted
                    ### Initial Conditions - not sensitive
                    a = np.array([[1,1,1,1]]).T
                    epsilon = 1
                    
                    # Resetting the error variables used to terminate the while loop
                    old_epsilon = 10
                    d_epsilon = 1
                    
                    # Set initial population distributions
                    S, Sv, E, I, R, In = initial_group_populations(N_vec,is_vacc,N_vec_vacc)
                    
                    count = 0
                    
                    # Loop fitting epsilon and transmission rates until converged
                    while d_epsilon > tol and (count < max_count):
                        # Update old epsilon
                        old_epsilon = epsilon
                        
                        
                        # Construct non-parametric matrix from census to compared to assortative matrix
                        C_unscaled = np.zeros([4,4])
                        for i in range(np.shape(eth_data)[0]):
                            C_unscaled += (eth_data[i:(i+1)].T@eth_data[i:(i+1)])*(a@a.T)/(eth_data[i:(i+1)]@a)
            
                        C_prime = C_unscaled
                        
                        # Optimise epsilon value
                        epsilon_soln = scipy.optimize.minimize_scalar(epsilon_variance,
                                             bounds=[0,1], args=(a,N_vector_eth, C_prime)
                                             )
                        epsilon = epsilon_soln.x

                        # optimise transmission rates
                        sol = scipy.optimize.least_squares(model_simulation_assortative,
                                                           a.flatten(),
                                                           bounds = (0.1,2), args=(CAR, np.concatenate([S,Sv,E,I,R,In]), N_vec, time, sigma,gamma,epsilon))
                        a = np.array([sol.x]).T
                        
                        # Calculate difference in epsilon values
                        d_epsilon = abs(epsilon - old_epsilon)
                        count += 1
                    
                    if count == max_count:
                        print('WARNING!'+'-'*81)
                        print('Fit terminated after {count} iterations because fit hit max iterations')
                        print('-'*90)
                    
                    save_epsilon_a_vector(epsilon,a,CAR,is_vacc,is_prop= False,is_non_parametric=False,
                                              is_SA1=is_SA1,is_statsnz=is_statsnz,is_old = False)
                    print(print_statement)
                    print()
                    print(f'Transmission rates {a}')
                    print(f'epsilon: {epsilon}')
                    print(f'error: {d_epsilon}  (want < {tol})')
                    print()


