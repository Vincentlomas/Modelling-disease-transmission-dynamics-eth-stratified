# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 15:17:34 2025

@author: Vincent Lomas

Module code
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import scipy


def epsilon_a_vector(CAR,is_SA1,is_statsnz,is_vacc,is_prop,is_non_parametric,is_old = False):
    '''Returns a tuple of a float, corresponding to the value of epsilon,
    and a 2D numpy array of shape (4,1), coresponding to the transmission rates'''
    
    if is_old:
        path = "transmission_rates/old_rates/"
    else:
        path = "transmission_rates/"
    
    if not os.path.exists(path): # if path doesn't exist, stop
        print('Path DNE')
    elif is_prop and is_non_parametric:
        print('Cannot be both proportionate mixing and non-parametric')
    else:
        if is_old:
            load_str0 = f'{path}old_'
        else:
            load_str0= path
        if is_vacc:
            load_str1 = load_str0 + 'vaccine_'
        else:
            load_str1 = load_str0
        if is_prop:
            load_str2 = load_str1 + f'prop_mix_CAR{CAR}.npy'
            loaded_array = np.load(load_str2)
            epsilon1 = loaded_array[0]
            a1 = np.array([loaded_array[1:]]).T
            return epsilon1, a1
        elif is_non_parametric:
            load_str2 = load_str1 + 'non-parametric_'
        else:
            load_str2 = load_str1 + 'assortative_'
        if is_statsnz:
            load_str3 = load_str2 + 'total_'
            if is_SA1:
                load_str4 = load_str3 + 'SA1_'
            else:
                load_str4 = load_str3 + 'SA2_'
        else:
            load_str4 = load_str2 + 'prioritised_SA2_'
        load_str5 = load_str4 + f"CAR{CAR}.npy"
        loaded_array = np.load(load_str5)
        epsilon1 = loaded_array[0]
        a1 = np.array([loaded_array[1:]]).T
        return epsilon1, a1


def filename_suffix(CAR,is_SA1,is_statsnz,is_non_parametric,is_prop,counterfactual):
    '''returns the suffix for filenames depending on information passed in'''
    if is_prop:
        suffix = 'prop_mix'
    else:
        if is_statsnz:
            suffix = 'total_measure'
        else:
            suffix = 'priority_measure'
    if not CAR is None:
        suffix = suffix + f'_CAR{CAR}'
    
    if is_non_parametric:
        suffix = 'non_parametric_' + suffix
    
    if not is_prop:
        if is_SA1 is None:
            suffix = suffix
        elif is_SA1:
            suffix += '_SA1'
        else:
            suffix += '_SA2'
        
    if counterfactual == 0:
        suffix += '_equal_max_vaccines'
    elif counterfactual == 1:
        suffix += '_no_vaccines'
    elif counterfactual == 2:
        suffix += '_SEIRS'
    elif counterfactual == 3:
        suffix += "_no_assortativity"
    elif counterfactual == 4:
        suffix += "_avged_transmission_rate"
    elif counterfactual == 5:
        suffix += "_avged_vaccination_rate"
        
    return suffix



def initial_group_populations(pop_vec,is_vacc,pop_vec_vacc=np.array([]),initial_exposed=0.0001,initial_infected=0,initial_recovered=0):
    '''takes as input:
        pop_vec: a vector of total populations of different ethnicities
        pop_vec_vacc: a vector of vaccinated populations of different ethnicities
        initial_exposed: the proportion of people initially exposed to the disease
        initial_infected: the proportion of people initially infected with the disease
        initial_recovered: the proportion of people initially recovered from the disease
    Returns initial SEIR group populations as a 5 length tuple of the form
    (susceptible, susceptible_vaccinated, exposed, infected, reovered)'''
    # Convert any 2D arrays to 1D
    pop_vec = pop_vec.flatten()
    pop_vec_vacc = pop_vec_vacc.flatten()
    
    initial_s = 1 - initial_exposed - initial_infected - initial_recovered
    
    if is_vacc:
        Sv = pop_vec_vacc
        S = initial_s*pop_vec - Sv
        E = initial_exposed*pop_vec
        I = initial_infected*pop_vec
        R = initial_recovered*pop_vec
    else:
        S = initial_s*pop_vec
        Sv = np.zeros(4)
        E = initial_exposed*pop_vec
        I = initial_infected*pop_vec
        R = initial_recovered*pop_vec
    In = np.zeros(4)
    return S,Sv,E,I,R, In



def load_beta(N_vec,CAR,is_vacc,is_prop,is_non_parametric,
                          is_SA1=None,is_statsnz=None, eth_data = None, is_old = False):
    
    epsilon,a = epsilon_a_vector(CAR,is_SA1,is_statsnz,is_vacc,is_prop,is_non_parametric,is_old)
    if is_non_parametric:
        b = non_parametric_per_capita_transmission_matrix(a, N_vec, eth_data = None)
    else:
        b = (1-epsilon)*(a @ a.T)/(a.T @ N_vec) + epsilon * np.diag(a[:,0] / N_vec[:,0])
    return b



def load_SEIR_results(CAR,is_vacc,is_prop,is_non_parametric,
                          counterfactual, is_SA1=None,is_statsnz=None,is_old = False):
    '''Returns the SEIR model run associate with the input in the order of:
        ts: a 1d array of times
        S: a 2d array of shape (4,len(ts)) of susceptible unvaccinated people
        Sv: a 2d array of shape (4,len(ts)) of susceptible vaccinated people
        E: a 2d array of shape (4,len(ts)) of exposed people
        I: a 2d array of shape (4,len(ts)) of infectious people
        R: a 2d array of shape (4,len(ts)) of recovered people
        In: a 2d array of shape (4,len(ts)) of the running total of infections'''
    
    if is_old:
        path = "SEIR_model_results/old_rates/"
    else:
        path = "SEIR_model_results/"
    
    if not os.path.exists(path): # if path doesn't exist, stop
        print('Path DNE')
    if is_prop and is_non_parametric:
        print('Cannot be both proportionate mixing and non-parametric')
    elif (is_SA1 == None or is_statsnz == None) and not is_prop:
        print('Define is_SA1 and is_statsnz when not loading proportionate mixing')
    else:
        if is_old:
            load_str0 = f'{path}old_SEIR_'
        else:
            load_str0= f'{path}SEIR_'
        if is_vacc:
            load_str1 = load_str0 + 'vaccine_'
        else:
            load_str1 = load_str0
        if is_prop:
            load_str2 = load_str1 + f'prop_mix_CAR{CAR}.npy'
            loaded_array = np.load(load_str2)
        else:
            if is_non_parametric:
                load_str2 = load_str1 + 'non-parametric_'
            else:
                load_str2 = load_str1 + 'assortative_'
            if is_statsnz:
                load_str3 = load_str2 + 'total_'
                if is_SA1:
                    load_str4 = load_str3 + 'SA1_'
                else:
                    load_str4 = load_str3 + 'SA2_'
            else:
                load_str4 = load_str2 + 'prioritised_SA2_'
            if counterfactual == -1:
                load_str5 = load_str4 + f"CAR{CAR}.npy"
            else:
                load_str5 = load_str4 + f"CAR{CAR}_counterfactual{counterfactual}.npy"
            loaded_array = np.load(load_str5)
        ts = loaded_array[0,:]
        S = loaded_array[1:5,:]
        Sv= loaded_array[5:9,:]
        E = loaded_array[9:13,:]
        I = loaded_array[13:17,:]
        R = loaded_array[17:21,:]
        In= loaded_array[21:,:]
        return ts, S, Sv, E, I, R, In



def model_populations(is_vacc):
    '''Returns a pair of vectors that contain the populations used for the model
    These are 4D arrays and returned in the order of: total ethnic populations,
    total vaccinated ethnic populations,
    Order of ethnicities is Maori, Pacific, Asian, European/Other'''
    
    mao_N = 888840
    pac_N = 359480
    asi_N = 820580
    oth_N = 3048470
    N_vec = np.array([[mao_N,pac_N,asi_N,oth_N]]).T

    ### vaccination populations
    # From data
    mao_N_vacc = 514930
    pac_N_vacc = 274413
    asi_N_vacc = 629841
    oth_N_vacc = 2558007
    N_vec_vacc = np.array([[mao_N_vacc,pac_N_vacc,asi_N_vacc,oth_N_vacc]]).T
    
    return N_vec, N_vec_vacc


def non_parametric_per_capita_transmission_matrix(a, N_vec, eth_data = None):
    '''Takes as input:
        a: the contact vector
        SA1: a boolean that is true if statistical area 1 data is needed and 
        false if SA2 data is needed
        is_statsnz: a boolean that is true for total data and flase for priority
        
        outputs: a per capita transmission matrix consturnced from proportionate
        mixing within statistical areas. martix outputted is 4 by 4.'''
    
    # only looking at SA2 prioritised, hence is_SA1 = False, is_statsnz=False
    if eth_data is None:
        eth_data = SA_import(is_SA1=False, is_statsnz=False)
    
    # Draw populations from vectors
    N_vector_eth = np.sum(eth_data,axis=0,dtype='int64')
    
    if np.all(N_vec.flatten()==N_vector_eth): # Check populations match
    
        # Caculating transmissions in SAs
        C_unscaled = np.zeros([4,4])
        for i in range(np.shape(eth_data)[0]):
            if not sum(eth_data[i]) == 0:
                C_unscaled += (eth_data[i:(i+1)].T@eth_data[i:(i+1)])*(a@a.T)/(eth_data[i:(i+1)]@a)
    
        # Some potential for normalising matrix if needed - We removed scaling from Ma et al.    
        C_prime = C_unscaled
        
        return (C_prime/(N_vector_eth)) /N_vec
    else: # Print error
        print('WARNING! '+81*'-')
        print("The population of the data set does not match what you input")
        print(90*'-')



### Proportionate mixing
def initial_reproduction_number(N_vec, N_vec_vacc,CAR,is_SA1,is_statsnz,is_vacc,is_prop,is_non_parametric, gamma):
    
        # Grab transmission rates and assortativity values
        epsilon, a = epsilon_a_vector(CAR,is_SA1,is_statsnz,is_vacc,is_prop,is_non_parametric,is_old = False)
        
        if is_non_parametric: # Handle non-parametric transmission matrix
            if is_statsnz:
                eth_data = SA_import(is_SA1, is_statsnz)
            else:
                eth_data = None
            beta = non_parametric_per_capita_transmission_matrix(a, N_vec, eth_data)

        else: # Handle proportionate and assortative transmission matrix
            beta = (1-epsilon)*(a @ a.T)/(a.T @ N_vec) + epsilon * np.diag(a[:,0] / N_vec[:,0])
        
        if is_vacc:
            beta_pop_matrix = beta*(N_vec - 0.3*N_vec_vacc)/gamma
        else:
            beta_pop_matrix = beta*(N_vec)/gamma
            
        
        F_Vinv = np.zeros([8,8])
        F_Vinv[0:4,0:4] = beta_pop_matrix
        F_Vinv[0:4,4:] = beta_pop_matrix
        
        eigs = np.linalg.eigvals(F_Vinv)
        
        return max(abs(eigs))
    


def return_attack_rates(CAR):
    if CAR == '4060':
        attack_rate_mao = 23.911 *100/ 40
        attack_rate_pac = 30.528 *100/ 40
        attack_rate_asi = 17.077 *100/ 60
        attack_rate_oth = 18.564 *100/ 60
    else:
        attack_rate_mao = 23.911 *100/ CAR
        attack_rate_pac = 30.528 *100/ CAR
        attack_rate_asi = 17.077 *100/ CAR
        attack_rate_oth = 18.564 *100/ CAR
    return np.array([attack_rate_mao,attack_rate_pac,attack_rate_asi,attack_rate_oth])



def save_epsilon_a_vector(epsilon1,a1,CAR,is_vacc,is_prop,is_non_parametric,
                          is_SA1=None,is_statsnz=None,is_old = False):
    '''saves a 5 length array from an input of an epsilon value and a transmission
    vector'''
    
    save_vector = np.zeros(5)
    save_vector[0] = epsilon1
    save_vector[1:] = a1.flatten()
    
    if is_old:
        path = "transmission_rates/old_rates/"
    else:
        path = "transmission_rates/"
    
    if not os.path.exists(path): # if path doesn't exist, create file and warn
        print('Path did not exist, created')
        os.makedirs(path)
    if is_prop and is_non_parametric:
        print('Cannot be both proportionate mixing and non-parametric')
    elif (is_SA1 == None or is_statsnz == None) and not is_prop:
        print('Define is_SA1 and is_statsnz when not saving proportionate mixing')
    else:
        if is_old:
            save_str0 = f'{path}old_'
        else:
            save_str0= path
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
            save_str5 = save_str4 + f"CAR{CAR}.npy"
            np.save(save_str5, save_vector)



def save_image(name_str,CAR,is_statsnz, is_SA1, is_prop, is_non_parametric, is_vacc, counterfactual,is_tight = True,file_format="png",dpi_val = 300):
    """Saves images to file location
    name_str is a string of the name of the file (without file extension)
    defaults to empty str, if empty string saves as new untitled image
    is_prop is a boolean, True if proportionate mixing, False otherwise, defaults
    to False
    CAR is case assertainment rate in percent, 40, 50, or 60, defaults to 50
    file_format is a string and is the extention type of the file, defaults to pdf
    dpi_val is an integer and is the dpi of the save image, defaults to 300"""
    
    # Make sure axis titles look nice
    if is_tight:
        plt.tight_layout()
    
    # Add vaccine to path if is_vacc
    if is_vacc:
        vacc_str = "Vaccine_analysis/"
        name_str = "vaccine_" + name_str
    else:
        vacc_str = ""
    
    if CAR in [40,50,60,'4060',None]: # Check Case assertainment rate is valid value
        # if CAR value, put in that folder
        if CAR is None:
            CAR_folder = ''
        else:
            CAR_folder = f'CAR{CAR}/'
        # Check if proportionate mixing or not
        if is_prop:
            # Define file path for image location
            file_path = f'../Images/{vacc_str}Proportionate_mixing/{CAR_folder}'
        elif is_non_parametric:
            # Define file path for image location
            file_path = f'../Images/{vacc_str}Non-parametric_mixing/{CAR_folder}'
        else:
            # Define file path for image location
            file_path = f'../Images/{vacc_str}Assortative_mixing/{CAR_folder}'
        
        if not counterfactual == -1:
            if counterfactual == 0:
                counter_fact_str = 'Equal_max_vaccines'
            elif counterfactual == 1:
                counter_fact_str = 'No_vaccines'
            elif counterfactual == 2:
                counter_fact_str = 'SEIRS_model'
            elif counterfactual == 3:
                counter_fact_str = 'No_assortativity'
            elif counterfactual == 4:
                counter_fact_str = 'Average_transmission_rates'
            elif counterfactual == 5:
                counter_fact_str = 'Averaged_vaccination_rates'
            file_path = f'../Images/Counterfactual_scenarios/{counter_fact_str}/'
        # If no name passed, save as untitled image
        if name_str == "":
            count=1
            name_str = f"untitled{count}_{filename_suffix(CAR,is_SA1,is_statsnz,is_non_parametric,is_prop, counterfactual)}"
            # Loop over folder to check for untitled names until a unique one is found
            while os.path.isfile(f'{file_path}{name_str}.{file_format}'):
                count += 1
                name_str = f"untitled{count}"
        else: # Just append suffix otherwise
            name_str = f'{name_str}_{filename_suffix(CAR,is_SA1,is_statsnz,is_non_parametric,is_prop, counterfactual)}'
        # If folder doesn't exist create it
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        # Save the image and print name and location
        plt.savefig(f'{file_path}{name_str}.{file_format}', dpi=dpi_val)
        print(f"Saved as {name_str}.{file_format} at {file_path}")
    else:
        print("CAR not 40, 50, 60, '4060', or None")
    print()



def SA_import(is_SA1, is_statsnz=True, file_location = '../Data', return_sa_ids = False):
    '''This function takes in two parameters: is_SA1 and is_statsnz
    is_SA1 has no default value and is a boolean, if true statistical area 1 
    data is imported, if false statistical area 2
    is_statsnz defaults to True and is a boolean, if true imports non-prioritised
    census data from statsnz if false imports Te Whatu Ora prioritised census
    data'''
    if is_statsnz:
        if is_SA1: # Filenames are defined for my system
            filename = f'{file_location}/SA1_ethnicity_data_2018.csv'
        else:
            filename = f'{file_location}/SA2_ethnicity_data_2018.csv'
            
        # Shape is (valid rows) by (4 ethnicities)
        # column names: Area_Code	European	Māori	Pacific Peoples	Asian	Middle Eastern / Latin American / African	Other Ethnicity	New Zealander(19)	Other Ethnicity nec(19)	Total stated	Not Elsewhere Included	 Total
        # New Zealander and Other Ethnicity nec are sub categories of other
        
        # import data
        raw_data = np.genfromtxt(filename, delimiter=',', encoding="utf8") 
        
        # Remove rows with invalid or confidential data in any position
        raw_data = raw_data[~np.ma.fix_invalid(raw_data).mask.any(axis=1)]
        
        # Create matrix with columns we care about, columns are: Māori	Pacific_Peoples	Asian	European/other
        c_data = np.zeros([np.shape(raw_data)[0],4],dtype="int64")
    
        # input relevant data
        c_data[:,0] = raw_data[:,2] # Maori
        c_data[:,1] = raw_data[:,3] # Pacific
        c_data[:,2] = raw_data[:,4] # Asian
        c_data[:,3] = raw_data[:,1] + raw_data[:,5] + raw_data[:,6] # Other = European + Middle_east + other
        
        # Grab SA ids
        SA_list = raw_data[:,0].astype('int32').tolist()
    else:
        if is_SA1:
            c_data = None # No SA1 data exists for the Te Whatu Ora data
        else:
            filename = f'{file_location}/SA2(Projections_as_at_30_Jun)_2022.csv'
            # columns are: sa2id	eth.mpao	popcount
            # import data
            raw_data = np.genfromtxt(filename,skip_header=1,dtype="i8,U8,U8,U7,U8,U8,i4,i8,U8,U8,U8", delimiter=',', encoding="utf8")
            
            # create list to assign index to each SA
            SA_list=[]
            # set up data set
            c_data = np.zeros([len(raw_data)//3,4],dtype='int64')
            # iterate over imported data rows
            for ID, ignore1, ignore2, eth, ignore4, ignore5, ignore6, population, ignore8, ignore9, ignore10 in raw_data:
                if not ID in SA_list: # append to list if ID not yet covered
                    SA_list.append(ID)
                idx = SA_list.index(ID) # grab index
                # put into matrix depending on ethnicity type
                if eth == 'Maori':
                    c_data[idx,0] += population
                elif eth == 'Pacific':
                    c_data[idx,1] += population
                elif eth == 'Asian':
                    c_data[idx,2] += population
                elif eth == 'Other':
                    c_data[idx,3] += population
    # Remove zero rows
    SA_list = (np.array(SA_list)[np.any(c_data[0:len(SA_list),:] > 0, axis=1)]).tolist()
    c_data = c_data[np.any(c_data > 0, axis=1)]
    # Output
    if return_sa_ids:
        return SA_list, c_data
    else:
        return c_data



def SEIR_model(t,SEIR, beta, sigma, gamma, imm_decay_rate=0):
    '''Takes as input a time (t), which is unused within the function, and
    a length 20 1D array the first 4 values coresspond to the susceptible 
    unvaccinated group for the 4 ethnic groups, next four to the suceptible
    vaccinated groups, and so on for E, I, and R as well
    The expected change in a timespan of 1 if conditions remain the same  for
    each value in the SEIR vector (i.e. the gradient) is calculated using the
    SEIR model
    The gradient of each value in the SEIR vector is returned resulting in a
    length 20 1D array'''
    S=SEIR[0:4] 
    Sv=SEIR[4:8]
    E=SEIR[8:12] 
    I=SEIR[12:16] 
    R=SEIR[16:20]
    In=SEIR[20:]
    ### Equations
    S_to_E = (beta.T *S).T @ I
    Sv_to_S = Sv * imm_decay_rate
    Sv_to_E = (beta.T *Sv).T @ I * (1-0.3) # Vaccination effectiveness
    E_to_I = sigma * E
    I_to_R = gamma * I
    R_to_S = R * imm_decay_rate # Make people susceptible again
    
    # Calculate migration gradients
    dS =  Sv_to_S + R_to_S - S_to_E
    dSv=- Sv_to_S - Sv_to_E
    dE =   S_to_E + Sv_to_E - E_to_I
    dI =   E_to_I - I_to_R
    dR =   I_to_R - R_to_S
    dIn =  E_to_I
    return np.concatenate([dS,dSv,dE,dI,dR,dIn])

