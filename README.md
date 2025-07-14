# Modelling-the-interaction-between-ethnicity-and-infectious-disease-transmission-dynamics
Repository to reproduce results for paper

# How to run code
## Update all results
In the SEIR_models folder there is a python file called "SEIR_model_update_all.py" 
In this file there are various parameters

time is the number of days the simulation is run for
sigma and gamma are model parameters that control the rate at which people move between the E, I, and R groups
is_vacc determines if vaccination effects are in the model

is_generate_transmission_rates, is_generate_SEIR_results, is_generate_plots, is_save_generated_plots are boolean variable that will run sections of the code to fit transmission rates, run SEIR models, generate plots, and save those plots respectively

## Changing datasets
If data sets are changed then they can be added to the SA_import function in the thesis_modules.py file and called upon when needed

## Statistical analysis of Case/Statistical area datasets
The timeseries_analysis.py file generates plots and stats about case data 

The SA_statistics_calculation.py file generates statistics about the statistical areas

## Image location
Images will all be stored in the Images folder
If vaccines are used in the plot, they will be stored in the Vaccine_analysis folder, which is similar to the root folder

If a counterfactual scenario is being plotted, the images will be saved in the Counterfactual scenario folder (over vaccine)

If applicable they will be stored in the relevant mixing folder

Finally, all folders (save couterfactual) have CAR folders for 4060, 40, 50, and 60, hich is where the images will go if relevant
