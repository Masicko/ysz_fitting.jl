################## Typical basic usage ######################
#
# use this funtion to call simple simulation of the experiments
#   ysz_fitting.simple_run()
#
# parameters:
#   TC, pO2, (bias=0.0) >>> the experimental setting
#         note: experimental data are for
#                 TC \in    {700, 750, 800, 850}
#                 pO2 \in   {0, 20, 40, 60, 80, 100}
#                 bias \in  {-1.0, -0.5, 0.0, 0.5, 1.0} but for 700 is missing {-1.0, 1.0}
#         note: You can write an array instead of 1 value for this parameters (as for prms_values described below)
#         WARNING: pO2 < 1.0e-5 => pO2 = 1.0e-5 ... regularization treashold
#   simulations     >>> =["EIS"] or =["CV"] or =["EIS", "CV"] for perform both methods for all experimental settings
#                       or now also ["CAP"] means differential capacitance computed with R0 = 0
#                       =["CV(f)"]   
#                               \-------TODO !!!
#                       
#   pyplot          >>> 0 = plot nothing, 1 = finall plot of the experiment, 2 = plot details
#   prms_names      >>> in array or tuple of strings specify, which parametr of the model should be changed
#   prms_values     >>> spedify, to which value it will be changed
#         note: You can change any of the parameters in YSZParameters class by typing its name, e.g. x_frac = 0.50
#         note: You can do small parametrical study by typing an array of values instead of one value
#              For example:  prms_names = ["DGA", "DGR", "DGO"]
#                            prms_values = [0.4, 0.5, 0.3]                       # 1 set of parameters
#                            prms_values = [0.4, [0.4, 0.5, 0.6], 0.3]           # 3 sets of parameters
#                            prms_values = [0.4, collect(0.4 : 0.1 : 0.6), 0.3]  # 3 sets of parameters
#   show_experiment >>> if the experimental setting is in the database, it shows fitting error and (if pyplot > 1) plot experimental data
#   fig_size        >>> standard PyPlot tuple (x,y)
#
# side-notes:
#   - this is using LoMA model
#   - expA \in {0, 1} determines, if the kinetics of the reaction will be switched to pure exponential form (from LoMA)
#   - "L" inductance
#   - TC = (700, 750, 800, 850)°C  => use => DD = ( 1.277, 2.92, 5.35, 9.05)e-13
#       - I know it should influence also the other parameters like \nu, but this is the first approximation to match overall resistance of the system
#
#
# #### update: CAP stuff ####
# 
# ["CAP"] simulation is also implemented
#   - overpotential = 0 is the equilibrium voltage !!!
#   - parametr R0 is set to 0 (but it can be overwriten to non-zero by prms_values)
#   - there are no experimental data, you must use "use_experiment=false" switch! 
#   >>> ["CAP"] Analytical version  - using direct capacitance analytical solution for defined voltage checknodes in simulation_struct
#
#   >>> ["CAP-CV"]     - using voltammetry simulation with default voltrate = 0.000001
#                      - calculate the capacitance from a rescaled current ... C = I / (voltrate * direction_of_sweep)
#                           
#                         
#                       
##############################################################

includet("../src/ysz_fitting.jl")

# example 1
ysz_fitting.simple_run(TC=850, pO2=[20, 40], simulations=["EIS", "CV"], pyplot=1, 
prms_names=["betaA", "betaR", "betaO", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[0.5, 0.5, 0.5,       0, 0, 0,       21.0, 20.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13])

# example 2 (this is equivalent to example 1)
ysz_fitting.simple_run(TC=850, pO2=[20, 40], simulations=["EIS", "CV"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      21.0, 20.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13])

# example 3
ysz_fitting.simple_run(TC=850, pO2=20, simulations=["EIS"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      collect(20.5 : 0.5 : 21.5), collect(19.5 : 0.5 : 20.5), 20.5,        0.7, 0.0, 0.0,     9.05e-13])

# example 4 (capacitance & not using the experimental data)
ysz_fitting.simple_run(TC=850, pO2=[20, 40], simulations=["CAP"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      21.0, 20.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13], use_experiment=false)

# example 5 (pO2 = 0 actually menas pO2 = 1.0e-5 as for all pO2 < 1.0e-5)
ysz_fitting.simple_run(TC=850, pO2=[0], simulations=["EIS", "CV", "CAP"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      21.0, 20.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13], use_experiment=false)

# example 6 (If for example R0 != 0 && K0 == 0, it provides higher capacitance)
ysz_fitting.simple_run(TC=850, pO2=[0], simulations=["CAP", "CAP-CV"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      21.0, 0.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13], use_experiment=false)

# example 7 - DRT analysis of two semi-circles
ysz_fitting.simple_run(TC=850, pO2=[80], simulations=["EIS"], pyplot=1,  
       prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu"     ], 
       prms_values=[1, 1, 1,       19.7, 19.8, 18.2,        0.0, -0.8, -0.0,     [90]*1.0e-14, [0.85]     ], 
       use_experiment=false);
       
# example 8 - DRT LoMA - two - circles .... near exp data
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_option="Nyq RC"), pyplot=1,  
              prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu"     ], 
              prms_values=[  23.6, 20.8, 20.2,            0.5, -0.4, 0.15, [90]*1.0e-14, [0.85]     ], 
              use_experiment=true);
              
# "OC" influence 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_option="Nyq RC"), pyplot=1,  
              prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "nus", "ms_par", "OC"     ], 
              prms_values=[  23.6, 20.8, 20.2,    0.5, collect(-0.4 : -0.05 : -0.4), 0.15,    [90]*1.0e-14, 0.85, 0.9, collect(0.05 : -0.01 : 0.05), 1 ./collect(3 : 0.4 : 8)      ], 
              use_experiment=true);

# 3 peakyyyyyy !!! :) - OVERVOLTAGE MODEL
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0), pyplot=1,  
                     prms_names=["e_fac", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu"     ], 
                     prms_values=[0.0,  collect(21.0 : 0.1 : 21), 19.7, 20.1,   0.8, -0.35, -0.3,    
                      [90]*1.0e-14, [0.85]     ], 
                     use_experiment=false);
# invariance test
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [0.001], 0.0, plot_option="RC"), pyplot=1,  
                     prms_names=[ "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu"     ], 
                     prms_values=[1,1,1,  collect(19 : 1.0 : 21), collect(19 : 1.0 : 21), collect(19 : 1.0 : 21),         0.0, -0.0, 0.0,     [90]*1.0e-14, [0.85]     ], 
                     use_experiment=false);
                     
# 3 peaky pro ysz_shn.jl model
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_option="DRT RC"), pyplot=1,  
                     prms_names=["kappaA", "kappaR", "kappaO", 
                                 "rA", "rR", "rO",     "rB", "rC",     
                                 "DGA", "DGR", "DGO",     
                                 "DD", "nu"     ], 
                     prms_values=(1.0, 0.0, 0.0,
                                  collect(21.65 : 0.05 : 21.65), collect(20.0 : -0.1 : 20.0), 19.7,     10, 10,       
                                  collect(-0.5 : 0.05 : -0.5), 0.0, 0.05,    
                                  [90]*1.0e-14, [0.85]     ), 
                     use_experiment=false);

# 4 peaky !!!!!!!
novy = 1; ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_option="Nyq Rf RC"), pyplot=1,  
                     prms_names=["kappaA", "kappaR", "kappaO", 
                                 "rA", "rR", "rO",     "rB", "rC",     
                                 "DGA", "DGR", "DGO",     
                                 "DD", "nu", "sites_Om0", "sites_Om1"     ], 
                     prms_values=(1.0, 0.0, 0.0,
                                  collect(21.65 : 0.05 : 21.65), collect(20.0 : -0.1 : 20.0), 19.7,     collect(18 : 0.05 : 20), 10,       
                                  collect(-0.5 : 0.05 : -0.5), 0.0, 0.05,    
                                  [90]*1.0e-14, [0.85], 2.0, collect(0.07 : 0.01 : 0.07)     ), 
                     use_experiment=false);
