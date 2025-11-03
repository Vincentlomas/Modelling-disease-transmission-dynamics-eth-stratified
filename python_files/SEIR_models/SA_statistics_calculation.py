# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 10:57:16 2025

@author: Vincent Lomas

Calculation of mean, median, SD, and range for statistical areas
"""

import pandas as pd
import numpy as np
from thesis_modules import SA_import

def print_stats(is_SA1, is_statsnz=True):
    '''Print in latex table format the total, mean, meadian, standard deviation,
    range, and inner quartile range'''
    ### Import data
    eth_import_data = SA_import(is_SA1, is_statsnz)
    eth_data = np.zeros([np.shape(eth_import_data)[0],5])
    eth_data[:,:4] = eth_import_data
    eth_data[:,4] = np.sum(eth_data,axis=1)
    
    eth_stats = pd.DataFrame(eth_data).describe()
    eth_median = np.median(eth_data, axis=0)
    
    print(f'M\=aori & {np.sum(eth_data[:,0]):.0f} & {eth_stats.iloc[1,0]:.1f} & {eth_median[0]} & {eth_stats.iloc[2,0]:.1f} & {eth_stats.iloc[3,0]:.0f} - {eth_stats.iloc[-1,0]:.0f} & {eth_stats.iloc[4,0]:.0f} - {eth_stats.iloc[6,0]:.0f}\\\\')
    print(f'Pacific Peoples & {np.sum(eth_data[:,1]):.0f} & {eth_stats.iloc[1,1]:.1f} & {eth_median[1]} & {eth_stats.iloc[2,1]:.1f} & {eth_stats.iloc[3,1]:.0f} - {eth_stats.iloc[-1,1]:.0f} & {eth_stats.iloc[4,1]:.0f} - {eth_stats.iloc[6,1]:.0f}\\\\')
    print(f'Asians & {np.sum(eth_data[:,2]):.0f} & {eth_stats.iloc[1,2]:.1f} & {eth_median[2]} & {eth_stats.iloc[2,2]:.1f} & {eth_stats.iloc[3,2]:.0f} - {eth_stats.iloc[-1,2]:.0f} & {eth_stats.iloc[4,2]:.0f} - {eth_stats.iloc[6,2]:.0f}\\\\')
    print(f'Europeans/Others & {np.sum(eth_data[:,3]):.0f} & {eth_stats.iloc[1,3]:.1f} & {eth_median[3]} & {eth_stats.iloc[2,3]:.1f} & {eth_stats.iloc[3,3]:.0f} - {eth_stats.iloc[-1,3]:.0f} & {eth_stats.iloc[4,3]:.0f} - {eth_stats.iloc[6,3]:.0f}\\\\')
    print(f'All ethnicities & {np.sum(eth_data[:,4]):.0f} & {eth_stats.iloc[1,4]:.1f} & {eth_median[4]} & {eth_stats.iloc[2,4]:.1f} & {eth_stats.iloc[3,4]:.0f} - {eth_stats.iloc[-1,4]:.0f} & {eth_stats.iloc[4,4]:.0f} - {eth_stats.iloc[6,4]:.0f}\\\\')

print_stats(True, is_statsnz=True)

# sa_ids_pri, eth_data_pri_import = SA_import(is_SA1=False, is_statsnz=False, return_sa_ids = True)
# sa_ids_tot, eth_data_tot_import = SA_import(is_SA1=False, is_statsnz=True, return_sa_ids = True)

# eth_data_pri = np.zeros([np.shape(eth_data_pri_import)[0],5])
# eth_data_pri[:,:4] = eth_data_pri_import
# eth_data_pri[:,4] = np.sum(eth_data_pri,axis=1)

# eth_data_tot = np.zeros([np.shape(eth_data_tot_import)[0],5])
# eth_data_tot[:,:4] = eth_data_tot_import
# eth_data_tot[:,4] = np.sum(eth_data_tot,axis=1)

# # Find unqies SA ids for each dataset
# unique_pri_ids = []
# for sa_id in sa_ids_pri:
#     if not sa_id in sa_ids_tot:
#         unique_pri_ids.append(sa_id)

# unique_tot_ids = []
# for sa_id in sa_ids_tot:
#     if not sa_id in sa_ids_pri:
#         unique_tot_ids.append(sa_id)



# sa_ids_pri_bigger = []
# sa_ids_pri_bigger_count = np.zeros(5)
# for i, sa_id in enumerate(sa_ids_pri):
#     if not sa_id in unique_pri_ids:
#         idx = sa_ids_tot.index(sa_id) # grab index
#         pri_greater = eth_data_pri[i,:] > eth_data_tot[idx,:]
#         sa_ids_pri_bigger_count += pri_greater
#         for j in np.array([0,1,2,3,4])[pri_greater]:
#             sa_ids_pri_bigger.append((j,sa_id))

# sa_ids_pri_bigger_amount_increase = np.zeros(5)
# for idx, sa_id in sa_ids_pri_bigger:
#     idx_tot = sa_ids_tot.index(sa_id)
#     idx_pri = sa_ids_pri.index(sa_id)
#     sa_ids_pri_bigger_amount_increase[idx] += eth_data_pri[idx_pri,idx] - eth_data_tot[idx_tot,idx]
# sa_ids_pri_bigger_amount_increase = sa_ids_pri_bigger_amount_increase/sa_ids_pri_bigger_count