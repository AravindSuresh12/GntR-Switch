This is a systems biology ML model which predicts responses for gene regulation in a genetic circuit. This mechanistic model predicts the transcriptional and translational dynamics present in the circuit.

WORK UNDER PROGRESS...


You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):


$ git pull https://github.com/AravindSuresh12/GntR-Full-Switch


or you can clone it from :

$ git clone https://github.com/AravindSuresh12/GntR-Full-Switch 

User is required to install Julia. This can be done at: https://julialang.org/. Necessary packages will be installed when user runs the file Include.jl from the terminal the first time. 


The science behind this work can be found with a similar work: https://www.biorxiv.org/content/10.1101/2023.01.10.523462v1.abstract

Files and their functions

|File|Description

| --- | --- |
	
|parameter_estimation_W_splined.jl|Solves the model equations for the ensemble of parameters sets for the test case, 10mM gluconate. Saves solutions in the poets_ensemble_W_test directory.|
| --- | --- |

Driver.jl- runs the simulated solution of ODEs present in the model.

Updated_Driver.jl- Gives a fit of the model's predicted solution versus experimental data obtained.

calculate_ensembles.jl- calculates the ensemble of ODE solutions for various concentrations of gluconate. The first 100 simulations of 2mM and 10mM are given in this repo as a reference.

Include.jl- run this file first to start the model. Installs all necessary packages if needed. 


