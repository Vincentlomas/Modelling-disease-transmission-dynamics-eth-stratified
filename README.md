# Modelling-the-interaction-between-ethnicity-and-infectious-disease-transmission-dynamics
Repository to reproduce results for paper

# How to run code
## Packages
Before running the code various module may need to be downloaded, these are: cycler, scipy, numpy, matplotlib, and seaborn

To do this you will need to install packages using [pip](https://packaging.python.org/en/latest/tutorials/installing-packages/) or equivelent.

## Update all results
In the SEIR_models folder there is a python file called "SEIR_model_update_all.py" when this file is run it will call on all other files except for Timeseries_analysis.py and SA_statistics_calculation.py. This file will update all result if specified using various parameters:

time is the number of days the simulation is run for
sigma and gamma are model parameters that control the rate at which people move between the E, I, and R groups
is_vacc determines if vaccination effects are in the model

is_generate_transmission_rates, is_generate_SEIR_results, is_generate_plots, is_save_generated_plots are boolean variables that will run sections of the code to fit transmission rates, run SEIR models, generate plots, and save those plots, respectively.

## Changing datasets
If data sets are changed then they can be added to the SA_import function in the thesis_modules.py file and called upon when needed. Some code (such as the SA_import module will need to have the path/name of the dataset added)

## Statistical analysis of Case/Statistical area datasets
The Timeseries_analysis.py file generates plots and stats about case data and needs to be run sperately from the other files (as it is unrelated to the SEIR model)

The SA_statistics_calculation.py file generates statistics about the statistical areas and formats them into a latex table for use in the paper

## Image location
Images will all be stored in the Images folder and code will generate the folders if needed
If vaccines are used in the plot, they will be stored in the Vaccine_analysis folder, which is similar to the root folder

If a counterfactual scenario is being plotted, the images will be saved in the Counterfactual scenario folder (instead of vaccine)
Note that this is where the quantification analysis will be present ()

If applicable they will be stored in the relevant mixing folder

Finally, all folders (save couterfactual) have CAR folders for 4060, 40, 50, and 60, hich is where the images will go if relevant


One such example of an image has been put as a text file in the location that an assortative mixing vaccine plot with a CAR of 40% would go
