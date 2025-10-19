# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:40:48 2024

@author: Vincent Lomas

A Python script to do statistical analysis and plots of Covid Timeseries data
"""

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large',
         'axes.titleweight': 'bold',
         'figure.titleweight': 'bold'}
pylab.rcParams.update(params)

colours=['#12436D', '#28A197', '#801650', '#F46A25']

### Key_word is admission, case, or death to get respective information
keyword = "case"

# Read in data
cdata = pd.read_csv(f"../Data/{keyword}s_by_eth/covid19_data.csv")
# Change dates to datetime data type
cdata[f'{keyword}_date'] = pd.to_datetime(cdata[f'{keyword}_date'])
cdata[f'{keyword}_week_ending'] = pd.to_datetime(cdata[f'{keyword}_week_ending'])

# ### Moving averages --- ALREADY IN DATASET
# cdata['MA_cases'] = cdata.apply(
#     lambda x: cdata[(cdata['ethnicity_mpao'] == x['ethnicity_mpao']) &
#                     (cdata['case_date'] <= x['case_date']) &
#                  (cdata['case_date'] > x['case_date'] - pd.Timedelta('7 days'))]['cases'].sum()/7,axis=1)

# Boolean vectors to seperate out ethnicities
mao_vec = cdata['ethnicity_mpao'] == "Maori"
pac_vec = cdata['ethnicity_mpao'] == "Pacific peoples"
asi_vec = cdata['ethnicity_mpao'] == "Asian"
oth_vec = cdata['ethnicity_mpao'] == "European or Other"

# Boolean vector to isolate the first omicron wave
omi1_vec = np.logical_and(cdata[f'{keyword}_date'] >= pd.to_datetime("2022-02-01"),
                    cdata[f'{keyword}_date'] <= pd.to_datetime("2022-05-16"))

# omi1_vec = np.logical_and(cdata[f'{keyword}_date'] >= pd.to_datetime("2022-05-16"),
#                     cdata[f'{keyword}_date'] <= pd.to_datetime("2022-10-16"))

# Boolean vector for full 2022-23
# omi1_vec = np.logical_and(cdata[f'{keyword}_date'] >= pd.to_datetime("2022-01-01"),
#                     cdata[f'{keyword}_date'] <= pd.to_datetime("2023-12-31"))     

# # Boolean vector for full COVID span up to 2023
# omi1_vec = np.logical_and(cdata[f'{keyword}_date'] >= pd.to_datetime("2020-02-01"),
#                     cdata[f'{keyword}_date'] <= pd.to_datetime("2023-12-31"))     


# Omicron wave 1 and population boolean vectors
mao_omi1_vec = np.logical_and(omi1_vec,mao_vec)
asi_omi1_vec = np.logical_and(omi1_vec,asi_vec)
pac_omi1_vec = np.logical_and(omi1_vec,pac_vec)
oth_omi1_vec = np.logical_and(omi1_vec,oth_vec)


mao_percent_omi1 = np.sum(cdata[f'{keyword}s_per_100k'].values[mao_omi1_vec]) / 1000
asi_percent_omi1 = np.sum(cdata[f'{keyword}s_per_100k'].values[asi_omi1_vec]) / 1000
pac_percent_omi1 = np.sum(cdata[f'{keyword}s_per_100k'].values[pac_omi1_vec]) / 1000
oth_percent_omi1 = np.sum(cdata[f'{keyword}s_per_100k'].values[oth_omi1_vec]) / 1000

x1 = cdata[f'{keyword}_date'].values[mao_omi1_vec]
y1 = cdata[f'{keyword}s_mean7'].values[mao_omi1_vec]
z1 = cdata[f'{keyword}s_per_100k'].values[mao_omi1_vec]
x2 = cdata[f'{keyword}_date'].values[pac_omi1_vec]
y2 = cdata[f'{keyword}s_mean7'].values[pac_omi1_vec]
z2 = cdata[f'{keyword}s_per_100k'].values[pac_omi1_vec]
x3 = cdata[f'{keyword}_date'].values[asi_omi1_vec]
y3 = cdata[f'{keyword}s_mean7'].values[asi_omi1_vec]
z3 = cdata[f'{keyword}s_per_100k'].values[asi_omi1_vec]
x4 = cdata[f'{keyword}_date'].values[oth_omi1_vec]
y4 = cdata[f'{keyword}s_mean7'].values[oth_omi1_vec]
z4 = cdata[f'{keyword}s_per_100k'].values[oth_omi1_vec]


cases1 = cdata[f'{keyword}s'].values[mao_omi1_vec]
cases2 = cdata[f'{keyword}s'].values[pac_omi1_vec]
cases3 = cdata[f'{keyword}s'].values[asi_omi1_vec]
cases4 = cdata[f'{keyword}s'].values[oth_omi1_vec]

mao_N =  802030 # From dataset cases/cases_per_100k
pac_N =  390970
asi_N =  834055
oth_N = 3167397

plt.rcParams.update({'font.size': 14})
fig, axs = plt.subplots(1,2,figsize=(6.8*2,4.8))

axs[0].plot(x1, y1, label = "Maori",color=colours[0])
axs[0].plot(x2, y2, label = "Pacific Peoples",color=colours[1])
axs[0].plot(x3, y3, label = "Asian",color=colours[2])
axs[0].plot(x4, y4, label = "European/other",color=colours[3])
axs[0].set_xlabel("Date")
axs[0].set_ylabel(f"daily {keyword}s MA")
axs[0].tick_params(labelrotation=-15,labelsize=10)
axs[0].set_title('a) Confirmed cases plot')
axs[0].legend(prop={'size': 10})
# plt.savefig(f"../Images/Raw_data_{keyword}s_7dayMA_COVID_start-2023.pdf",dpi=300)
# plt.show()

axs[1].plot(x1, y1 / mao_N * 100000, label = "Maori",color=colours[0])
axs[1].plot(x2, y2 / pac_N * 100000, label = "Pacific Peoples",color=colours[1])
axs[1].plot(x3, y3 / asi_N * 100000, label = "Asian",color=colours[2])
axs[1].plot(x4, y4 / oth_N * 100000, label = "European/other",color=colours[3])
axs[1].set_xlabel("Date")
axs[1].set_ylabel(f"daily {keyword}s MA (per 100k)")
axs[1].tick_params(labelrotation=-15,labelsize=10)
axs[1].legend(prop={'size': 10})
axs[1].set_title('b) Per capita confirmed cases plot')
plt.tight_layout()
plt.savefig(f"../Images/Raw_data_{keyword}s_feb1-may16.png",dpi=300)
plt.show()
plt.rcParams.update({'font.size': 10})

# plt.xlabel("Ethnicity")
# plt.ylabel("Percent of ethnicity")
# plt.title("Covid cases in first Omicron wave (Feb 1st - Apr 16th 2022)")
# plt.savefig("../Images/Percent_cases_by_eth_omi1.png",dpi=300)
# plt.show()

# fig = plt.figure()
# # creating the bar plot
# plt.bar(['Maori','Pacific peoples','Asian', 'European/other'],
#         [mao_percent_omi1, pac_percent_omi1, asi_percent_omi1, oth_percent_omi1],
#         color ='maroon', width = 0.8)

# # plt.xlabel("Ethnicity")
# plt.ylabel("Confirmed cases per ethnicity (%)")
# plt.tight_layout()
# # plt.savefig("../Images/barplot_attack_rate_data.png",dpi=300)
# plt.show()
