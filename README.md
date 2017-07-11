# SABS_first_project

--------------------------------------------------------------------------------------------------------
																									   
   "Inversion based on simultaneous observations of voltage and calcium concentration in iPSC-CMs"    
																									   
                           Joanneke E Jansen (joanneke.jansen@maths.ox.ac.uk)						   
                                             April-June 2017										   
																									   
--------------------------------------------------------------------------------------------------------

This repository contains the directories 'Latex', 'Code' and 'Various'. It also contains my final report and two powerpoint presentations: final.pdf, final_short.pptx and final_long.pptx.

1) 'Latex'
Contains my report. Old versions of the report to be found in the folder 'Old versions'.

2) 'Code'
Contains four directories: 'Tissue', 'Single Cell', 'Extra' and 'cbcbeat'. The 'Single Cell' directory contains all code for single cell simulations, the 'Tissue folder' contains all code for the monodomain simulations. We took a 12mm strip as our domain for the Tissue simulations.

	2a) 'Tissue'
	* create_initial_conditions_12mm_strip.py can be used to create an h5 file: initial_conditions_12mm_strip.h5. The model is run for 800s.

	* thereafter, create_observations_12mm_strip.py must be used to create synthetic observations: observed_v_A_B_C_D_12mm_strip.h5, observed_cai_A_B_C_D_12mm_strip.h5, and observed_times_A_B_C_D_12mm_strip.txt. Here, A,B,C and D stand for the value of the factor we multiply g_Na, g_Kr, g_K1 and g_CaL with respectively. Those values should be given as arguments and should lie in the range [0.5-1.5].

	* alternatively, create_noisy_observations_12mm_strip.py can be used to create noisy observations. The noise percentages can be given as arguments. The output files contain the noise percentages in their names.

	* inverse_12mm_strip.py solves the inverse problem. The levels of noise can be passed as arguments. For example, when the code is run with the command
	"python inverse_12mm_strip -gnaf 1.0 -gkrf 1.0 -gk1f 1.0 -gcalf 1.0", synthetic observations with g_Na, g_Kr, g_K1 and g_CaL equal to their default values are used. 

	2b) 'Single Cell'
	* similar to the 'Tissue' code.

	2c) 'Extra'
	* files to calculate the values of the functionals efficiently, while varying various parameters at a time.

	* files to save our results in a format that can be viewed in Paraview.

	2d) 'cbcbeat'
	* the 'cbcbeat' folder contains the Paci2013 cell models we use. 

3) 'Various'
Contains various other files. The folders 'Matlab' and 'Paraview' contain the data and scripts to generate the figures in the report.