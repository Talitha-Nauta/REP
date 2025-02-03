# Simulation for delay attacks on cyber pysical systems

## Files
main.py defines the systems and global parameters, it runs the other files using exec(optimise_{}.py).  
optimise_bin_chi2.py does the optimisation using both the binominal and the chi2 test.
optimise_bin.py does the optimisation using the binominal test.
optimise_chi2.py does the optimisation using the chi2 test.

The results are stored in the folders results_{bin/chi2/bin_chi2}/{system} and for each test configuration there are
two files named a{a_max}_p{p_miss}_alpha{significance_level}_w{window_length}. One csv file and one txt file. 
The csv file contains all information needed to make to plots, the txt files contains extra information such as the sequences.

main_short.py does the same as main.py but does only perform the test on system = furuta_zero, a_max = 3, p_miss = 0.25, window_length = 20

plot.tex creates a pdf with all the figures in the paper.
plot_short.tex creates a pdf with a figure for the short version.
The makefile can be used to render the pdf, comment the corresponding line for normal/short.
It is possible to create the plots even when the code is not finished with the longer sequences.

## Versions
python          3.12.2
gurobipy        11.0.1
numpy           1.26.4
scipy           1.12.0

## Install Gurobi and obtain a license
To install gurobi and obtain a license follow the instruction on 
https://www.gurobi.com/features/academic-named-user-license/  

## Setting up a virtual environment
If your python installation support venv to create virtual environments, the following commands would make it possible to create a virtual environment and run the code inside said virtual environment on a linux. This requires that you have a gurobi license, otherwise only parts of the execution will work.   

% cd REP
comment: navigate to the folder where you saved the reproducibility package REP

% python3 -m venv .
comment: this creates the folder bin and initializes the virtual environment

% source ./bin/activate
comment: this activates the virtual environment,
        your terminal shows (REP) at the beginning from now on
        
(REP) % python3 -m pip install gurobipy numpy scipy
comment: this installs the packages in the virtual environment

(REP) % python3 main_short.py
comment: expect this command to take 2 hours


