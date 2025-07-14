# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 11:17:02 2025

@author: Vincent Lomas

Main file that will update all results based on input datasets
"""

from thesis_modules import *
from non_parametric_matrix_parameter_optimisation import non_parametric_optimisation
from assortative_mixing_determination import assortative_optimisation
from proportionate_mixing_optimisation import proportionate_optimisation
from SEIR_model_code import run_SEIR_model
from plot_code import *

# Model parameters
time = 365 # length of time to fit and run parameters
sigma = 1/3 # Rate of disease development
gamma = 0.25 # Recovery rate
is_vacc = True # determines if vaccination effects are in the model

is_generate_transmission_rates = False
is_generate_SEIR_results = False
is_generate_plots = True
is_save_generated_plots = True

### Populations being worked with
N_vec, N_vec_vacc = model_populations(is_vacc)
# # Old Populations - not for use, but for interpretation of old results
# mao_N =  802030 # From dataset cases/cases_per_100k
# pac_N =  390970
# asi_N =  834055
# oth_N = 3167397
# tot_N = mao_N + asi_N + pac_N + oth_N
# N_vec = np.array([[mao_N],
#          [pac_N],
#          [asi_N],
#          [oth_N]])

if is_generate_transmission_rates:
    print('Proportionate optimisation')
    # Fit proportionate mixing transmission rates
    proportionate_optimisation(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma)
    
    print('Assortative optimisation')
    # Fit assortative mixing transmission rates
    assortative_optimisation(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma)
    
    print('non-parametric optimisation')
    # Fit non-parametric mixing transmission rates
    non_parametric_optimisation(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma)

if is_generate_SEIR_results:
    
    # Proportionate mixing
    run_SEIR_model(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma, is_prop=True, is_non_parametric=False)
    
    # Assortative mixing
    for is_SA1, is_statsnz in [(True,True),(False,True),(True,True)]:
        run_SEIR_model(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma, is_prop=False,
                       is_non_parametric=False, is_SA1 = is_SA1, is_statsnz = is_statsnz)
        
    # Non-parametric mixing
    run_SEIR_model(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma, is_prop=False,
                   is_non_parametric=True, is_SA1 = False, is_statsnz = False)
    
    
    ### counterfactual
    # -1 is no scenario, 0 is equal vaccine, 1 is no vaccines, 2 is SEIRS model
    # 3 is no assortativity
    # 4 is setting all contact rates to weighted average contact rate
    # 5 is setting all vaccination rates to weighted average vaccine rate
    for counterfactual in [0,3,4,5]:
        for is_SA1, is_statsnz in [(True,True),(False,True),(True,True)]:
            run_SEIR_model(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma, is_prop=False,
                           is_non_parametric=False, is_SA1 = is_SA1, is_statsnz = is_statsnz, counter_factual_scen = counterfactual)
    for counterfactual in [1,2]:
        run_SEIR_model(N_vec, N_vec_vacc, is_vacc, time, sigma, gamma, is_prop=False,
                       is_non_parametric=True, is_SA1 = False, is_statsnz = False, counter_factual_scen = counterfactual)
        
if is_generate_plots:
    CAR = 50
    is_SA1 = False
    is_statsnz = False
    
    ### SEIR plot of each ethnicity
    is_prop = True
    is_non_parametric = False
    counterfactual = -1
    SEIR_plot(N_vec,CAR,is_vacc,is_prop,is_non_parametric,is_SA1,is_statsnz,counterfactual)
    if is_save_generated_plots:
        save_image('per_capita_SEIR_plots',CAR,is_statsnz, is_SA1,is_prop, is_non_parametric, is_vacc,counterfactual,file_format = 'png')
    
    
    ### Comparrison of contact rates and initial reproduction number
    is_prop = True
    is_non_parametric = False
    counterfactual = -1
    is_SA1 = None
    is_statsnz = None
    contact_vs_reproduction_number(N_vec, N_vec_vacc,is_SA1,is_statsnz,
                                       is_vacc,is_prop,is_non_parametric, gamma=gamma)
    if is_save_generated_plots:
        save_image('transmission_rate_and_reproduction_number',CAR=None,is_statsnz=is_statsnz,
                   is_SA1=is_SA1,is_prop=is_prop, is_non_parametric=is_non_parametric,
                   is_vacc=is_vacc,counterfactual=counterfactual)
    
    ### Comparrison of contact rates
    contact_vs_contact_rates(N_vec, N_vec_vacc,is_SA1s = [False,False],
                             is_statsnzs=[False,False],
                             is_vacc=is_vacc,
                             is_props = [False,False],
                             is_non_parametrics = [False,True],
                             titles=["Assortative","Non-parametric"]) # transmission rates will follow
    if is_save_generated_plots:
        save_image('transmission_rate_and_reproduction_number',CAR=50,is_vacc=is_vacc,is_prop=False,is_non_parametric=True,is_SA1=is_SA1,is_statsnz=is_SA1,counterfactual=-1,file_format = 'png')
    
    ### Heat plots
    # Comparing all three mixing methods
    beta_prop = load_beta(N_vec,CAR,is_vacc,is_prop=True,is_non_parametric=False)
    beta_assort = load_beta(N_vec,CAR,is_vacc,is_prop=False,is_non_parametric=False,is_SA1=is_SA1,is_statsnz=is_statsnz)
    beta_non_para = load_beta(N_vec,CAR,is_vacc,is_prop=False,is_non_parametric=True,is_SA1=is_SA1,is_statsnz=is_SA1)
    heat_map((beta_prop*N_vec,beta_assort*N_vec,beta_non_para*N_vec,beta_prop,beta_assort,beta_non_para),['Proportionate','Assortative','Non-parametric','Per capita proportionate','Per capita assortative','Per capita non-parametric'],
              scaling=[0,0,0,7,7,7],create_fig=True,labels=["Māo","Pac","Asi","Oth"],max_comp='rowwise')
    if is_save_generated_plots:
        save_image('heatplot_transmission_proportionate_vs_assortative_vs',CAR=50,is_vacc=True,is_prop=False,is_non_parametric=True,is_SA1=is_SA1,is_statsnz=is_SA1,counterfactual=-1,file_format = 'png')
    
    # Comparing total measure fits of assortative mixing
    is_statsnz = True
    beta_SA1 = load_beta(N_vec,CAR,is_vacc,is_prop=False,is_non_parametric=False,is_SA1=True,is_statsnz=is_statsnz)
    beta_SA2 = load_beta(N_vec,CAR,is_vacc,is_prop=False,is_non_parametric=False,is_SA1=False,is_statsnz=is_statsnz)
    heat_map((beta_SA1*N_vec,beta_SA2*N_vec,beta_SA1,beta_SA2),['SA1','SA2','Per capita SA1','Per capita SA2'],
              scaling=[0,0,7,7],create_fig=True,labels=["Māo","Pac","Asi","Oth"],max_comp='rowwise')
    if is_save_generated_plots:
        save_image('heatplot_transmission_SA1_vs_SA2',CAR=50,is_vacc=True,is_prop=False,is_non_parametric=False,is_SA1=None,is_statsnz=is_SA1,counterfactual=-1,file_format = 'png')
    
    
    ### Quantification analysis
    CAR = 50
    for is_SA1, is_statsnz in [(True,True),(False,True),(False,False)]:
        quantification_plot(N_vec, CAR, is_SA1, is_statsnz)
        save_image('quantification_analysis', CAR, is_statsnz, is_SA1, is_prop=False,is_non_parametric=False, is_vacc=True, counterfactual=3,file_format='png')
        
        
        
    ### Epsilon variation
    epsilon_variation_plot(CAR,N_vec,is_vacc,N_vec_vacc, round(1.5*time), sigma,gamma,is_save_generated_plots)
    
    
