# -*- coding: utf-8 -*-
"""
Created on Thu May 15 11:31:19 2025

@author: Vincent Lomas
"""

from thesis_modules import *
import numpy as np
import os
import scipy


def save_SEIR_results(SEIR,times,CAR,is_vacc,is_prop,is_non_parametric,
                          counterfactual, is_SA1=None,is_statsnz=None,is_old = False):
    '''saves the results from an SEIR model run'''
    
    save_vector = np.zeros([len(SEIR)+1,len(times)])
    save_vector[0,:] = times
    save_vector[1:,:] = SEIR
    
    if is_old:
        path = "SEIR_model_results/old_results/"
    else:
        path = "SEIR_model_results/"
    
    if not os.path.exists(path): # if path doesn't exist, create file and warn
        print('Path did not exist, created')
        os.makedirs(path)
    if is_prop and is_non_parametric:
        print('Cannot be both proportionate mixing and non-parametric')
    elif (is_SA1 == None or is_statsnz == None) and not is_prop:
        print('Define is_SA1 and is_statsnz when not saving proportionate mixing')
    else:
        if is_old:
            save_str0 = f'{path}old_SEIR_'
        else:
            save_str0= f'{path}SEIR_'
        if is_vacc:
            save_str1 = save_str0 + 'vaccine_'
        else:
            save_str1 = save_str0
        if is_prop:
            save_str2 = save_str1 + f'prop_mix_CAR{CAR}.npy'
            np.save(save_str2, save_vector)
        else:
            if is_non_parametric:
                save_str2 = save_str1 + 'non-parametric_'
            else:
                save_str2 = save_str1 + 'assortative_'
            if is_statsnz:
                save_str3 = save_str2 + 'total_'
                if is_SA1:
                    save_str4 = save_str3 + 'SA1_'
                else:
                    save_str4 = save_str3 + 'SA2_'
            else:
                save_str4 = save_str2 + 'prioritised_SA2_'
            if counterfactual == -1:
                save_str5 = save_str4 + f"CAR{CAR}.npy"
            else:
                save_str5 = save_str4 + f"CAR{CAR}_counterfactual{counterfactual}.npy"
            np.save(save_str5, save_vector)



def run_SEIR_model(N_vec, N_vec_vacc, N_vec_vacc_boosted, is_vacc, time, sigma, gamma, is_prop, is_non_parametric, is_SA1 = None, is_statsnz = None, counter_factual_scen = -1, CARs = ['4060',40,50,60]):
    
    if (counter_factual_scen in [0,1,5]) and not is_vacc:
        counter_factual_scen = -1
        print(f'Warning is_SA1: {is_SA1}, is_statsnz: {is_statsnz} -- changed counterfactual scenario to -1')
    elif (is_prop or is_non_parametric) and (counter_factual_scen == 3):
        counter_factual_scen = -1
        print(f'Warning is_SA1: {is_SA1}, is_statsnz: {is_statsnz} -- changed counterfactual scenario to -1')
    elif counter_factual_scen == 0:
        # Manually setting equal to highest rate (european)
        N_vec_vacc = N_vec * (N_vec_vacc[-1,0] / N_vec[-1,0])
        N_vec_vacc_boosted = N_vec * (N_vec_vacc_boosted[-1,0] / N_vec[-1,0])
    elif counter_factual_scen == 1:
        # Setting to zero
        N_vec_vacc = np.zeros([4,1])
        N_vec_vacc_boosted = np.zeros([4,1])
    elif counter_factual_scen == 5:
        N_vec_vacc = sum(N_vec_vacc)/sum(N_vec) * N_vec
        N_vec_vacc_boosted = sum(N_vec_vacc_boosted)/sum(N_vec) * N_vec

    ### iterate over all case ascertainment rates
    for CAR in CARs:
    # Isolate cases if needed, if all cases to be run set to if True:
        if (counter_factual_scen == -1) or (is_vacc and CAR ==50):
            epsilon, a = epsilon_a_vector(CAR,is_SA1,is_statsnz,is_vacc,is_prop,is_non_parametric,is_old = False)
            
            if counter_factual_scen == 3:
                epsilon = 0
            elif counter_factual_scen ==4:
                a = np.full((4,1),sum(a*N_vec)/sum(N_vec))
            
            # Grab attack rates
            attack_rates = return_attack_rates(CAR)
                
            ### Assign beta, the transmission matrix
            if is_non_parametric:
                # Construct matrix using SA
                eth_data = SA_import(is_SA1, is_statsnz)
                beta = non_parametric_per_capita_transmission_matrix(a, N_vec, eth_data)
            else:
                # Construct using equation
                beta = (1-epsilon)*(a @ a.T)/(a.T @ N_vec) + epsilon * np.diag(a[:,0] / N_vec[:,0])
             
            ### initial group population distributions
            # 0.01% initial expoure for each group
            S, Sv, Svb, E, Ev, Evb, I, Iv, Ivb, R, In = initial_group_populations(N_vec,is_vacc,N_vec_vacc, N_vec_vacc_boosted)
                
            SEIR_0 = np.concatenate([S,Sv,Svb,E, Ev, Evb, I, Iv, Ivb,R,In])
            
            if counter_factual_scen == 2: # SEIRS model
                time_SEIRS = 4000 # length of time to run simulation (days)
                solution = scipy.integrate.solve_ivp(SEIR_model,[0,time_SEIRS],SEIR_0,args=(beta,sigma,gamma, 1/365),t_eval=np.arange(0,time+1))
            else:
                if counter_factual_scen == -1:
                    time = 365 # length of time to run simulation (days)
                else: # Run counter factuals for longer as some epidemics run longer than a year
                    time =365*2
                solution = scipy.integrate.solve_ivp(SEIR_model,[0,time],SEIR_0,args=(beta,sigma,gamma),t_eval=np.arange(0,time+1))
            
            
            # Define an acceptable error in attack rate for alerting if something is wrong
            tolerable_err = 0.001
            # Checking if the attack rates (in percentage units) are above the error threshold
            above_tol = abs(100*(N_vec.flatten() - solution.y[0:4,-1] - solution.y[4:8,-1] - solution.y[8:12,-1])/(N_vec.flatten()) - attack_rates.flatten()) > tolerable_err
            if np.any(above_tol) and counter_factual_scen == -1:
                print("-"*90)
                print(f"Warning CAR{CAR}, is_SA1: {is_SA1}, is_statsnz: {is_statsnz}, is above threshold")
                print("-"*90)
            
            if is_prop:
                save_SEIR_results(solution.y,solution.t,CAR,is_vacc,is_prop,is_non_parametric,
                                      counter_factual_scen, is_SA1=None,is_statsnz=None)
            else:
                save_SEIR_results(solution.y,solution.t,CAR,is_vacc,is_prop,is_non_parametric,
                                      counter_factual_scen, is_SA1,is_statsnz)
            
            if is_prop:
                print_statement = 'Prop mix: '
            elif is_non_parametric:
                print_statement = 'Non-para: '
            else:
                print_statement = 'Assort : '
            
            if counter_factual_scen == -1:
                print_statement_cf = ''
            else:
                print_statement_cf = f'counter factual: {counter_factual_scen}, '
            
            print(f"{print_statement}CAR{CAR}, is_SA1: {is_SA1}, is_statsnz: {is_statsnz}, {print_statement_cf}complete")
