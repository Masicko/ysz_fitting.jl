
# Temperature testing

new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [ 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.0, 0.0, 0.0, 0.0, 
                   0.0, 24.4626,            0.0, 0.33,
                   0.0, 22.63,              0.0, 0.0088,
                   0.0, 21.94,              0.0, 0.17,
                   0.0, 0.85,       0.0, 3.95,      0.0, 6.818*0.15)
               ,use_experiment=true);


# with step from 700 to 750               
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(700, [ 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 
                   0.2, 21.92,            -0.05, 0.055,
                   0.3, 21.21,             0.14, 0.048,
                   0.2, 21.809,             0.0, -0.04,
                   -0.1, 0.76,       20.0, 37.97,      4.0, 5.56*0.15)
               ,use_experiment=true);               
               
# fitted TC-all pO2-40-60 <><><><><><><><><><> ALL <><><>< HURAAAAA  ----> almost also CV fitted 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [40, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 
               0.22236513700278862, 21.92,               0.00465798980401979, 0.055, 
               0.33958803427894996, 21.21,              -0.06332811104733665, 0.048, 
               0.33062408382716696, 21.809,              0.15698240545620798, -0.04, 
               -0.253637378507218, 0.76,        9.999491268858524, 37.97,         0.21670880659050948, 0.834)
               ,use_experiment=true);
               
               
# fitted TC-all pO2-20-60 <><> ALL <><> HURAAAAA ------> EIS great !!! even 850 3-peaks, CV worse
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 
               0.2577645962065092, 21.89500309270786,       -0.02108143562365405, -0.003636121874040238, 
               0.37191997547659544, 21.249644272226337,      -0.04214106521785248, -0.005224773396108435, 
               0.3917042938923518, 21.56797364229525,          0.17238782294738939, -0.005598977167630661, 
               -0.1483966936929819, 0.76,         5.0497715462086745, 15.566557872263587,         0.20448043798844218, 0.6571331057177057)
               ,use_experiment=true);

# ~~~~~~~~~~~~~# temperature study show
                new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(collect(650 : 10 : 900), [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                        
                                      prms_values=(1, 0.27, 0.0, 0.0, 0.0, 
                                      0.2577645962065092, 21.89500309270786,       -0.02108143562365405, -0.003636121874040238, 
                                      0.37191997547659544, 21.249644272226337,      -0.04214106521785248, -0.005224773396108435, 
                                      0.3917042938923518, 21.56797364229525,          0.17238782294738939, -0.005598977167630661, 
                                      -0.1483966936929819, 0.76,         5.0497715462086745, 15.566557872263587,         0.20448043798844218, 0.6571331057177057)
                                      ,use_experiment=false);

# ~~~~~~~~~~~~~~# pO2 study show
                new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700], [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.27, 0.0, 0.0, 0.0, 
                      0.2577645962065092, 21.89500309270786,       -0.02108143562365405, -0.003636121874040238, 
                      0.37191997547659544, 21.249644272226337,      -0.04214106521785248, -0.005224773396108435, 
                      0.3917042938923518, 21.56797364229525,          0.17238782294738939, -0.005598977167630661, 
                      -0.1483966936929819, 0.76,         5.0497715462086745, 15.566557872263587,         0.20448043798844218, 0.6571331057177057)
                      ,use_experiment=false);
               
               
               
               
################# ############################ ############################## #######################
# EIS-CV fitting


[0.356, 0.304, 21.866, -0.0, 0.04, 0.378, 21.119, -0.05, 0.066, 0.258, 24.664, 0.212, -0.043, -0.175, 0.783,       7.922, 23.458, 0.408, 1.062]
[0.392, 0.476, 21.807, 0.075, 0.027, 0.283, 21.118, 0.055, 0.094, 0.009, 23.635, 0.091, -0.026, -0.293, 0.607, 2.958, 99.529, 0.629, 1.103]
[0.457, 0.371, 21.813, 0.039, 0.059, 0.3, 21.086, 0.022, 0.173, -0.007, 25.973, 0.115, -0.044, -0.139, 0.017, 2.128, 48.162, 1.105, 0.729] 
[0.082, 0.395, 21.597, 0.099, -0.095, 0.056, 22.582, -0.08, 0.042, 0.285, 21.243, -0.093, -0.043, -0.246, 0.018, 4.35, 1.723, 4.824, 1.476] 
        
###########################################        #############
###  ##                ######### # # 3  ###############################
# pO2 fittovani                           
#ff ["EIS", "CV"] 
#ff [1, 0.01]  -> [0.271,   0.236, 21.968, 0.009, 0.037, 0.349, 21.223, -0.062, 0.057, 0.44, 21.685, 0.188, -0.036, -0.254,   9.998, 30.35, 0.199, 0.856]
#ff [1, 0.1]   -> [0.296,   0.229, 21.965, 0.004, 0.007, 0.35, 21.228, -0.061, 0.078, 0.432, 21.668, 0.186, -0.035, -0.254,   9.97, 25.496, 0.218, 0.887]
#ff [1, 5]     -> [0.403,   0.239, 21.956, -0.004, 0.073, 0.278, 21.096, -0.056, 0.08, 0.731, 21.823, 0.151, -0.043, -0.254,  9.985, 51.946, 0.19, 0.797]
#ff [1, 10]    -> [0.36,    0.332, 21.916, 0.007, 0.011, 0.263, 21.148, -0.062, 0.084, 0.72, 21.848, 0.153, -0.047, -0.252,    9.999, 44.549, 0.13, 0.898]
#ff [1, 20]    -> [0.365,   0.375, 21.886, 0.0, 0.11, 0.282, 21.12, -0.059, 0.068, 0.609, 22.173, 0.162, -0.041, -0.264,      9.959, 19.805, 0.057, 1.095]
#ff [1, 100]   -> [0.341,   0.441, 21.683, 0.004, 0.073, 0.322, 21.102, -0.061, 0.044, 0.33, 23.26, 0.18, -0.041, -0.261,     9.732, 31.78, 0.161, 1.213]


#ff [1, 0.1]
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.296, 0.0, 0.0, 0.0, 
               0.229, 21.965, 0.004, 0.007, 0.35, 21.228, -0.061, 0.078, 0.432, 21.668, 0.186, -0.035, -0.254,
               
               0.76,
               
               9.97, 25.496, 0.218, 0.887)
               ,use_experiment=true); 

#ff [1, 5]
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.403, 0.0, 0.0, 0.0, 
               0.239, 21.956, -0.004, 0.073, 0.278, 21.096, -0.056, 0.08, 0.731, 21.823, 0.151, -0.043, -0.254,
               
               0.76,
               
               9.985, 51.946, 0.19, 0.797)
               ,use_experiment=true); 

#ff [1, 100]
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.341, 0.0, 0.0, 0.0, 
               0.441, 21.683, 0.004, 0.073, 0.322, 21.102, -0.061, 0.044, 0.33, 23.26, 0.18, -0.041, -0.261,
               
               0.76,
               
               9.732, 31.78, 0.161, 1.213)
               ,use_experiment=true);               
               
# GOOD both ... EIS two peaks ... but ok.. for pO2 20   ... BUT !!! simulation fails for higher pO2         
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.082, 0.0, 0.0, 0.0, 
               0.395, 21.597, 0.099, -0.095, 0.056, 22.582, -0.08, 0.042, 0.285, 21.243, -0.093, -0.043, -0.246, 0.018, 4.35, 1.723, 4.824, 1.476)
               ,use_experiment=true);                

               
               
               
               
               
               
               
               
               
               
               
               
               
#### ONLY EIS >>>>>>>>>>>> no CV fitting ###########################################               
#######################################################  COOL ###     but CV is bad         
# simulation for different pO2 20, 40
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 0.31448487403944203, 21.92349379172513, -0.016450793003275933, -0.076705719267254, 0.3203843257111922, 21.278403971316955, 0.013483504657377793, 0.0706847107563436, 0.16213347852298987, 21.48118352940071, 0.08111234036292242, -0.07871949753844011, -0.04999999999995275, 0.7168561989664564, -1.7230708878839274, 9.77546137619688, 0.34768745580475163, 0.6252830672992037)
               ,use_experiment=true);                

# test_pO2-40-60 ... only EIS   ----> CV good in OER and bad in ORR
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [40, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
              
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 0.26496545524910625, 22.00945665225656, 0.03502102061099427, 0.16812809974771498, 0.42687185198943034, 21.323795429611657, -0.10228956705169508, -0.049999999999903, 0.48660342436284626, 21.40125037246591, 0.19939370407133936, 0.027041695169748105, -0.24049176916593748, 0.7379769874108393, 9.993334722221602, 8.790844036559363, 0.24841207861994652, 0.6127277699241381)
               ,use_experiment=true);

               
# JEEEEEE ... pro pO2 20 je hezke EIS a dostatecne CV            
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 850], [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 0.357, 22.148, -0.038, 0.037, 0.256, 21.352, 0.11, 0.202, 0.435, 21.822, -0.089, -0.05, -0.05, 0.39, 2.123, 39.203, 2.912, 1.614)
               ,use_experiment=true); 
               
               
               
               
               
               
               
               
               
               
               
               
               
               
               
               
               

# long long fitting

new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.382, 0.0, 0.0, 0.0, 0.347, 21.731, 0.004, 0.041, 0.452, 21.136, -0.075, 0.086, 0.16, 25.13, 0.264, 0.009, -0.034, 0.027, -6.529, 64.263, 0.457, 1.148)
               ,use_experiment=true);   
               
               
# second hope ... ff[1, 0.5] -> CV bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.27, 0.0, 0.0, 0.0, 0.228, 22.133, -0.049, 0.13, 0.205, 21.337, 0.046, 0.081, 0.244, 21.745, -0.002, 0.061, -0.18, 0.42, 6.134, 16.788, 0.827, 0.776)
               ,use_experiment=true);
   
   
   
   
   

########### COOOOOOOOOOOOOOOOOOOOOOOOL ########   
# all_700, 850 pO2 60 --- ff[1 ,1]
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.394, 0.0, 0.0, 0.0, 0.291, 21.912, 0.035, -0.009, 0.208, 21.181, 0.017, 0.11, 0.57, 21.381, -0.006, -0.025, -0.149, 0.716, 2.085, 10.173, 0.251, 0.549)              
               ,use_experiment=true); 

# all - || - ... number two
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(1, 0.421, 0.0, 0.0, 0.0, 0.298, 21.959, -0.013, 0.141, 0.202, 21.203, 0.082, -0.05, 0.427, 21.353, -0.014, 0.019, -0.15, 0.745, 0.808, 9.136, 0.232, 0.579)              
               ,use_experiment=true);
               


# some temporary thing
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(0, 0.295, 0.0, 0.0, 0.0, 0.222, 21.836, 0.011, 0.138, 0.207, 24.047, -0.002, -0.075, 0.222, 21.74, 0.119, -0.05, -0.106, 0.703, 0.009, 28.533, 0.005, 1.853)               
               ,use_experiment=true);   

               
               
               
### TC 700-750-850, pO2-20-40-60       ---->>> unfinished because of not enough iterations !!!        
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(0, 0.463, 0.0, 0.0, 0.0, 0.519, 21.974, -0.09, 0.051, 0.387, 21.112, -0.074, -0.05, 0.25, 21.985, 0.016, -0.05, -0.078, 0.373, -0.024, 41.992, -0.011, 2.115)               
               ,use_experiment=true); 


### TC 700-750-850, pO2-20-40-60       ->>> finished ---- .. BUT !!! strange nu
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750, 800, 850], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
               prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                 
               prms_values=(0, 0.411, 0.0, 0.0, 0.0, 0.209, 22.366, -0.031, 0.174, 0.295, 20.997, -0.025, 0.227, 0.53, 21.859, 0.121, -0.19, 0.018, 0.387, 0.027, 57.057, -0.004, 1.925)
               ,use_experiment=true);

               
               
#####################
#####################   NEW      HOPE      ! !! ! ! ! ! ! ! !  comming fron diferential preliminary results

### standard SNEH configuarion
#  prms_lists=(0, 0.3, 0.0, 0.0, 0.0,
#                   collect(0.25 : 0.2 : 0.25), 21.8,            collect(-0.07 : 0.7 : 0.07), collect(-0.05 : 0.1 : 0.05),
#                   collect(0.25 : 0.2 : 0.25), 21.3,            collect(-0.07 : 0.7 : 0.07), collect(-0.05 : 0.1 : 0.05),
#                   collect(0.25 : 0.2 : 0.25), 21.5,            collect(-0.07 : 0.7 : 0.07), collect(-0.05 : 0.1 : 0.05),
#                   -0.15, 0.80,       5.0, 15,      0.2, 2)   #5.56*0.15)



new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=      (1, 0.32, 0, 0, 0,
                           0.2414, 21.9222,       -0.0618,     0.0558328, 
                           0.2818, 21.2101,       0.133,     0.0480154,
                           0.2607, 21.809,        -0.00615,    -0.0438488,
                           -0.133,  0.765282,     13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807) 
                      ,use_experiment=true);
               
               
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=      (1, 0.32, 0, 0, 0,
                           0.2414, 21.9222,       -0.0618,     0.0558328, 
                           0.2818, 21.2101,       0.133,     0.0480154,
                           0.2607, 21.809,        -0.00615,    -0.0438488,
                           -0.133,  0.765282,     13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807) 
                      ,use_experiment=true);               
               
###### RESULT #######
# ff EIS-CV = [1, 0.0] >>>>>>>> EIS very good
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.43096402436244524, 0, 0, 0, 0.288753339609282, 21.848922443832656, -0.027911417199004526, -0.004572027928900416, 0.34126793851080767, 21.125910388311805, -0.01811856204892673, 0.07188764966775679, 0.2812740178140478, 21.467392253748145, 0.13475333193229766, -0.05952919838333402, -0.11167034872750412, 0.6732907278249488, 0.0019313676794959437, 11.553528276411772, 0.1858817137725142, 0.5255525854475593)
                      ,use_experiment=true);

# ff EIS-CV = [1, 0.8] >>>>>>>> EIS very good
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.302, 0, 0, 0, 0.288, 21.868, -0.015, 0.033, 0.389, 21.199, -0.05, 0.049, 0.348, 21.754, 0.184, -0.022, -0.15, 0.661, 0.251, 33.803, 0.142, 0.848)
                      ,use_experiment=true);                    
                      
# ff EIS-CV = [1, 1.2] >>>>>>>> GOOD    BOTHT .... UA UA UA ... :) hura!  >>>>>>>> CORRECTED CO down
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.336, 0, 0, 0, 
                      0.332, 21.817,              -0.017, 0.024, 
                      0.397, 21.171,              -0.027, 0.042, 
                      0.937, 21.858,               0.195, -0.001, 
                      -0.15, 0.65,      9.996, 30.591,       0.713, 0.878)
                      ,use_experiment=true);
               
               
# ff EIS-CV = [1, inf] >>>>>>>> GOOD CV
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.2849959838598055, 0, 0, 0, 
                      0.45858880392144824, 21.874974975180162, 0.039069632415346334, 0.030162414682266617, 
                      0.3095420370564378, 21.07891961773628, -0.0008658579059998266, -0.04603628357043432, 
                      0.9911991334369225, 22.11768963100269, 0.10146097638438026, 0.004871994799156309, 
                      -0.11578447441283925, 0.9347985214936189, 14.715260842936846, 9.185710773669479, 3.7862431522855418, 1.8109281332760347)
                      ,use_experiment=true);                
               
               
               
               
               
               
               
               
############################# half-corrected ZETA #####################################################
##################################################################################################
# ff EIS-CV = [1, 1] 
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 0.275, 22.096, -0.035, 0.054, 0.351, 21.253, -0.05, 0.108, 0.336, 21.677, 0.163, -0.059, -0.144, 0.65, 2.615, 21.368, 0.202, 0.854)
                                                 ,use_experiment=true);
                                                 

# ff EIS-CV = [1, 2] 
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 0.31, 22.087, -0.021, 0.06, 0.337, 21.253, -0.017, 0.117, 0.797, 21.688, 0.162, -0.055, -0.15, 0.65, 22.577, 23.806, 0.868, 0.895)
                                                 ,use_experiment=true);
                                                 
                                               

# ff EIS-CV = [1, 3]
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 0.301, 22.125, -0.035, 0.123, 0.328, 21.245, -0.018, 0.127, 0.823, 21.695, 0.157, -0.06, -0.15, 0.65, 29.986, 26.939, 0.794, 0.923)                            
                                                 ,use_experiment=true);                                                 

# ff EIS-CV = [1, 4] ->>> good.... could be taken ... but A peak is weird ... :(
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 0.323, 22.054, -0.022, 0.101, 0.333, 21.243, -0.019, 0.104, 0.875, 21.707, 0.164, -0.051, -0.15, 0.65, 29.891, 32.197, 0.701, 0.977)                            
                                                 ,use_experiment=true); 
                  
                  
# ff EIS-CV = [1, 4 NONO] >> zeta = 0.35 !!!!!!!! NO NO NO NO 
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.35, 0, 0, 0, 0.339, 21.981, -0.035, 0.065, 0.394, 21.271, -0.05, 0.051, 1.141, 21.894, 0.206, -0.011, -0.104, 0.65, 28.266, 99.877, 0.571, 1.199)
                                                 ,use_experiment=true);                
             
# ff EIS-CV = [1, 5]
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 0.356, 22.04, -0.009, 0.115, 0.35, 21.227, -0.022, 0.104, 0.932, 21.721, 0.18, -0.055, -0.15, 0.65, 29.993, 73.423, 0.612, 1.107)
                                                 ,use_experiment=true);             
               
# ff EIS-CV = [1, 7]
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 0.364, 22.027, -0.009, 0.121, 0.351, 21.22, -0.024, 0.11, 0.944, 21.724, 0.184, -0.064, -0.147, 0.65, 15.52, 99.965, 0.554, 1.2)
                                                 ,use_experiment=true);
         
         
         
         
         
         
         
         
###############  useable single_TC guesses ### >> zeta = 0.39
# 700 
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700], pO2=[20, 60], data_set="OLD_MONO_100"          
          , pyplot=1,  
          prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp", 
              "A.r",      "A.DG",        
              "R.r",      "R.DG", 
              "O.r",      "O.DG",                      
              "nu",     "CO",     "COmm"],
          prms_values=(1, 0.39, 0, 0, 0, 22.2557, 0.196224, 21.2795, 0.148953, 21.7135, -0.0663272, 0.3, 18.4293, 0.845544)
          ,use_experiment=true);

# 750                                                                                                      
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[750], pO2=[20, 60], data_set="OLD_MONO_100"          
          , pyplot=1,  
          prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp", 
              "A.r",      "A.DG",        
              "R.r",      "R.DG", 
              "O.r",      "O.DG",                      
              "nu",     "CO",     "COmm"],
          prms_values=(1, 0.39, 0, 0, 0, 22.2426, 0.000406918, 21.636, 0.0284407, 22.6145, 0.165328, 0.317062, 149.999, 1.7805) 
          ,use_experiment=true);                                                 
                                                 
               
               
# difference ->> (with corrected nu 0.3 -> 0.8)               
prms_diff = (0, 0, 0, 0, 0, -0.013100000000001444, -0.195817082, 0.3565000000000005, -0.1205123, 0.9009999999999998, 0.2316552, 0.01706200000000002, 131.56969999999998, 0.934956)
               
## guess in the temperature framework ##               
ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100"                
                , pyplot=1,  
                prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                prms_values=(1, 0.45, 0, 0, 0, 
                              -0.0131, 22.2557,          -0.195817082, 0.196224, 
                              0.3565, 21.2795,           -0.1205, 0.148953, 
                              0.9, 21.7135,           0.2316, -0.0663272,
                              0.017, 0.8,        20.00, 18.4293,         0.93, 0.845544)
                                                 ,use_experiment=true);               
               
               
               
               
               
               
################################################################################################################################
################################################################################################################################               
################    FINALLY CORRECT ! ! !!!  ! !! ! !   ########################################################################
# ff EIS-CV = [1, 0.0] >>>>>>>> EIS very good
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([700, 750], [20, 60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100"), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.43096402436244524, 0, 0, 0, 0.288753339609282, 21.848922443832656, -0.027911417199004526, -0.004572027928900416, 0.34126793851080767, 21.125910388311805, -0.01811856204892673, 0.07188764966775679, 0.2812740178140478, 21.467392253748145, 0.13475333193229766, -0.05952919838333402, -0.11167034872750412, 0.6732907278249488, 0.0019313676794959437, 11.553528276411772, 0.1858817137725142, 0.5255525854475593)
                      ,use_experiment=true);
               
               
               
# ff EIS-CV = [1, 1.8] >>>>>>> BOTH good :))))))))))))))))
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.4, 0, 0, 0, 0.409, 22.171, -0.05, 0.011, 0.166, 0.405, 22.58, -0.071, -0.048, 0.512, 0.389, 21.22, 0.103, 0.005, -0.296, -0.05, 0.689, 27.145, 16.145, 0.174, 0.684)
                      ,use_experiment=true);

# results for  ./temporary/CORRECT__TC-700-750-850_pO2-20-60__12h_with-S_ff1-18/   
res1 = (1, 0.4, 0, 0, 0, 0.409, 22.171, -0.05, 0.011, 0.166, 0.405, 22.58, -0.071, -0.048, 0.512, 0.389, 21.22, 0.103, 0.005, -0.296, -0.05, 0.689, 27.145, 16.145, 0.174, 0.684)  
res2 = (1, 0.4, 0, 0, 0, 0.246, 21.852, -0.03, 0.141, -0.06, 0.397, 22.167, -0.05, -0.027, 0.407, 0.608, 22.057, 0.12, 0.027, 0.044, -0.05, 0.651, 8.235, 65.056, 0.141, 0.883)  ########## ASI nejlpesi .. . . .  
res3 = (1, 0.4, 0, 0, 0, 0.242, 22.131, -0.05, 0.121, 0.102, 0.33, 21.822, -0.038, 0.032, 0.335, 0.385, 21.519, 0.102, -0.029, -0.217, -0.05, 0.77, 0.038, 64.336, 0.167, 0.793)   
res4 = (1, 0.4, 0, 0, 0, 0.276, 21.736, -0.013, 0.144, -0.153, 0.348, 21.666, -0.017, 0.024, 0.271, 0.765, 22.244, 0.116, 0.006, 0.117, -0.05, 0.651, 2.582, 56.543, 0.207, 0.887)

               
               
               
                 

######## TESTING res2 branch >>>>>>                 
# ff EIS-CV = [1, 0.5] >>>>>>> EIS good ... CV not bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.4, 0, 0, 0, 0.274, 21.867, -0.03, 0.073, 0, 0.265, 22.072, 0.028, 0.035, 0.4, 0.351, 21.592, 0.037, -0.036, 0, -0.05, 0.65, 18.995, 13.388, 0.346, 0.589)
                      ,use_experiment=true);                
                
# ff EIS-CV = [1, 1.9] >>>>>>> BOTH good :))))))))))))))))
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.4, 0, 0, 0, 
                        0.2825409348052233, 21.887472830668642,     -0.023307013734441284, 0.0980508181649717,    0, 
                        0.2508668451347317, 22.090559844352843,     0.03792340204410957, -0.004860699297832186,   0.4, 
                        0.3868217597699285, 21.870010619930834,     0.028391243693966506, 0.007005115247570771,   0, 
                        -0.04999999999980735, 0.650006103286414,    29.996, 44.0577,    0.2711, 0.7615)
                      ,use_experiment=true);         
               
 # ff EIS-CV = [1, 50] >>>>>>> CV hm...not bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.4, 0, 0, 0, 0.694, 21.988, -0.05, 0.094, 0, 0.109, 24.664, -0.05, -0.595, 0.4, 0.103, 23.884, -0.05, 0.799, 0, -0.05, 0.65, 10.169, 7.766, 1.608, 2.276)
                      ,use_experiment=true);              
               
 # ff EIS-CV = [1, 50] >>>>>>> CV hm...not bad . . . . INTERESTING !!!
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values= (1, 0.4, 0, 0, 0, 1.0185215132822192, 23.182674107082356, 0.13440543573183525, -0.39997828296758825, 0, 0.5761228751304889, 23.162955498506758, -0.04823389345853183, -0.5681798094492313, 0.4, 0.3312458201313922, 21.691782616050627, 0.011559093251970432, -0.051549137252482446, 0, -0.04999999999985964, 0.9219231046985951, 10.42800243707937, 97.32851834184535, 1.5785439673767248, 0.08274263700524723)
                      ,use_experiment=true);                
            
               
############# HHH ###########   UUU ############# RRR ################# AAAAA  ######  !!!!!!!! ###############              
########   HAND GUESS for log_10 (X.S) = (0, 0.2, 0)
# ff EIS-CV = [1, 1.8] >>>>>>> CV hm...not bad . . . . INTERESTING !!!
# new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
#                       pyplot=1,  
#                       prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
#                           "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
#                           "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
#                           "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
#                           "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
#                         
#                       prms_values= (1, 0.307, 0, 0, 0, 0.293, 21.917, -0.027, 0.082, 0.0, 0.339, 21.541, -0.007, 0.029, 0.2, 0.913, 21.968, 0.109, -0.011, 0.0, -0.05, 0.65, 19.545, 60.431, 0.416, 0.86)
#                       ,use_experiment=true);
               
               
               
               
################ pO2 = [20, 40, 60]
# ff EIS-CV = [1, 0] >>>>>>>  
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 40, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.3, 0, 0, 0, 0.278, 21.896, -0.005, -0.071, 0.0, 0.312, 21.616, -0.024, 0.072, 0.2, 0.246, 21.41, 0.105, -0.1, 0.0, -0.05, 0.65, 0.0, 8.54, 0.165, 0.552)
                      ,use_experiment=true, line_color_idx=2, evaluate_temp=true);
                      
# ff EIS-CV = [1, 0] >>> Free NU and COmm
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20,40, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_700", "nu_750", "nu_800", "nu_850",       "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values= (1, 0.3, 0, 0, 0, 0.246, 21.927, 0.005, -0.095, 0.0, 0.328, 21.589, -0.009, 0.059, 0.2, 0.18, 21.51, 0.093, -0.088, 0.0, 0.729, 0.572, 0.566, 0.202, -2.725, 13.798, 0.102, 0.668)
                      ,use_experiment=true, line_color_idx=1, evaluate_temp=true, save_dir="freeNU");
                      
##########################################################################
##########################################################################

                      
#### pO2 = [20, 60] #################                     
# ff EIS-CV = [1, 0] >>>>>>> EIS >>>  GOOD :)  >>>>>>>>>>>>>> HUSTY !!!! trefi to celkem pekne i dalsi atmosfery kysliku 40
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values= (1, 0.3, 0, 0, 0, 0.279, 21.883, 0.007, -0.061, 0.0, 0.298, 21.621, -0.029, 0.087, 0.2, 0.245, 21.434, 0.108, -0.11, 0.0, -0.05, 0.652, 0.0, 8.88, 0.158, 0.59)
                      ,use_experiment=true, line_color_idx=1, evaluate_temp=true);
                      
                      
# ff EIS-CV = [1, 1.8] >>>>>>> both >>>  both GOOD :)
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values= (1, 0.3, 0, 0, 0, 0.295, 21.915, -0.027, 0.082, 0.0, 0.34, 21.543, -0.006, 0.026, 0.2, 0.961, 21.991, 0.108, -0.008, 0.0, -0.05, 0.65, 8.466, 69.399, 0.417, 0.873)
                      ,use_experiment=true);               
              
              
# ff EIS-CV = [1, inf] >>>>>>>  CV  >>> so ... so ... 
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 800, 850], pO2=[20, 60], data_set="OLD_MONO_100", 
                      pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values= (1, 0.3, 0, 0, 0, 0.669, 21.989, -0.098, 0.111, 0.0, 0.515, 22.627, -0.14, -0.33, 0.2, 0.563, 21.82, 0.046, 0.013, 0.0, -0.05, 0.772, 1.588, 0.468, 1.488, 0.724)
                      ,use_experiment=true);              
                
               
               

               
               
               
               
               
               
               
               
               
               
#p1 = (1, 0.4, 0, 0, 0, 0.409, 22.171, -0.05, 0.011, 0.166, 0.405, 22.58, -0.071, -0.048, 0.512, 0.389, 21.22, 0.103, 0.005, -0.296, -0.05, 0.689, 27.145, 16.145, 0.174, 0.684)   
p1 = (1, 0.3, 0, 0, 0, 0.295, 21.915, -0.027, 0.082, 0.0, 0.34, 21.543, -0.006, 0.026, 0.2, 0.961, 21.991, 0.108, -0.008, 0.0, -0.05, 0.65, 8.466, 69.399, 0.417, 0.873)

p2 = (1, 0.4, 0, 0, 0, 
                        0.2825409348052233, 21.887472830668642,     -0.023307013734441284, 0.0980508181649717,    0, 
                        0.2508668451347317, 22.090559844352843,     0.03792340204410957, -0.004860699297832186,   0.4, 
                        0.3868217597699285, 21.870010619930834,     0.028391243693966506, 0.007005115247570771,   0, 
                        -0.04999999999980735, 0.650006103286414,    29.996, 44.0577,    0.2711, 0.7615)               
               
                                                             
                                                        
# E^LSM = zeta E^YSZ
# naive try with fixed zeta = 0.3, SR = 0.2 ... as fot the previous case
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                           physical_model_name="ysz_model_GAS_LoMA_generic", 
                             pyplot=1,
                             prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                 "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                 "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                 "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                 "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                               
                             prms_values=(1, 0.3, 0, 0, 0, 0.373, 22.14, -0.08, -0.1, 0.0, 0.301, 21.741, 0.083, 0.05, 0.2, 0.301, 21.739, -0.08, 0.101, 0.0, 0.804, 0.702, 0.6, 0.537, 10.038, 30.718, 0.716, 1.004) 
                             ,use_experiment=true);
                                                        
# hope ff 1:1.5 ... not the best, but 3 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                      
                                    prms_values=(1, 1.015, 0, 0, 0, 0.263, 21.618, 0.079, -0.081, -0.29, 0.32, 21.482, 0.125, 0.084, 0.317, 0.31, 22.043, 0.075, -0.082, 0.001, 0.814, 0.717, 0.6, 0.511, 10.127, 30.321, 0.641, 0.873)
                                    ,use_experiment=true);
                                                        
                                                
# hope ff 1:1.5 ... not the best, but 3 peaks >> better than the one above
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                      
                                    prms_values=(1, 0.942, 0, 0, 0, 0.323, 21.544, 0.063, -0.082, -0.326, 0.306, 21.527, 0.128, 0.087, 0.322, 0.35, 22.12, 0.056, -0.091, 0.006, 0.829, 0.7, 0.6, 0.496, 10.141, 32.648, 0.348, 0.698)
                                    ,use_experiment=true);                                                        
                                                        
###########################  FINALLLLLLLLLLLLLLLLLLLLYYYYYYYY 
# hope ff 1:0 ... great :))            >>> taken                          
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                      
                                    prms_values=(1, 0.88, 0, 0, 0, 
                                          0.217, 21.566,      0.073, -0.136,        -0.33, 
                                          0.277, 21.655,      0.025, 0.118,          0.33, 
                                          0.124, 21.291,      0.076, -0.133,         0.0, 
                                          0.75, 0.704, 0.521, 0.201,         -1.344, 6.574, -0.025, 0.395)                                    
                                    ,use_experiment=true, 
                                    save_files=false, save_dir="EIS", line_color_idx=1, evaluate_temp=true);      
                                    
# hope ff 1:0  >> CO(mm)_B nonneg                                    
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                      
                                    prms_values=(1, 0.88, 0, 0, 0,  
                                          0.264, 21.545,        0.045, 0.045,        -0.33, 
                                          0.306, 21.585,        0.012, 0.17,          0.33, 
                                          0.243, 21.291,        0.095, -0.1,           0.0, 
                                          0.785, 0.662, 0.47, 0.209,           0.163, 7.149, 0.089, 0.354)                                    
                                    ,use_experiment=true, 
                                    save_files=false, save_dir="EIS_nn", line_color_idx=1, evaluate_temp=true);                                    
                                                
                                                        
# hope ff 1:1.3 ... not the best, but 3 peaks >> GOOD >>> could be taken !!! :)
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                      
                                    prms_values=(1, 0.88, 0, 0, 0,
                                          0.338, 21.533,      0.067, -0.08,       -0.33, 
                                          0.333, 21.584,      0.144, 0.084,        0.33, 
                                          0.352, 21.983,      0.048, -0.087,       0.001, 
                                          0.819, 0.708, 0.6, 0.517,       9.989, 30.319,      0.486, 0.688)
                                    ,use_experiment=true, 
                                    save_files=false, save_dir="both", line_color_idx=2, evaluate_temp=true);

               

               
               
               
###### CONST CO(mm) ###########################################
# hope ff 1:0 ... great :))            >>> taken                          
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                      
                                    prms_values=(1, 0.88, 0, 0, 0, 0.257, 21.517, 0.043, -0.094, -0.33, 0.237, 21.732, -0.004, 0.202, 0.33, 0.18, 21.215, 0.085, -0.156, 0.0, 0.75, 0.679, 0.496, 0.2, 0.0, 4.116, 0.0, 0.38)                                 
                                    ,use_experiment=true,
                                    save_files=false, save_dir="EIS", line_color_idx=1, evaluate_temp=true);
               
               
               
# hope ff 1:1.3 ... 2 peaks ... but  >> GOOD >>> could be taken !!! :)
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                        
                                    prms_values=(1, 0.88, 0, 0, 0, 0.3, 21.6, 0.08, 0.11, -0.33, 0.3, 21.6, 0.103, 0.011, 0.33, 0.306, 22.307, 0.08, 0.0, 0.0, 0.8, 0.746, 0.646, 0.5, 0.0, 30.591, 0.0, 0.878)
                                    ,use_experiment=true, 
                                    save_files=false, save_dir="both", line_color_idx=2, evaluate_temp=true);               
# hope ff 1:1.3 ...    CV ok                
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=[700, 750, 850], pO2=[20, 60], data_set="OLD_MONO_100",
                                  physical_model_name="ysz_model_GAS_LoMA_generic", 
                                    pyplot=1,
                                    prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                                        "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                                        
                                    prms_values=(1, 0.88, 0, 0, 0, 0.54, 21.502, 0.013, 0.097, -0.33, 0.326, 22.234, 0.112, 0.113, 0.33, 0.426, 21.965, 0.093, -0.088, 0.0, 0.93, 0.65, 0.6, 0.26, 0.0, 30.591, 0.0, 0.878)
                                    ,use_experiment=true, 
                                    save_files=false, save_dir="CV", line_color_idx=3, evaluate_temp=true);               
               
               
               
               
               
               
               
               
               
               
               
               
               
        
        
        
        
        
        
        
        
        
        
        
        
        
               
               
           

                                              

### Conductivities #####
#
# OLD_MONO_100 -> [700, 750, 800, 850] = [ ???, ???,  3.58  , 5.75]
#
#
#
#

#######################################################################################################################################
################################################# separate ############################################################################
#######################################################################################################################################                                          
                                          
############################################################################
#################### EEE ###################################################
############################################################################

################# no_e_fac ###################
# EEE EIS CV ... relatively good fit for both CV and EIS
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=100, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                         prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(true, 1.0, 1.0, 1.0, 21.94051183287039, 21.553329968593776, 20.965273571151613, 0.09841737344413018, -0.091875316601428, 0.04652480321433385, 9.3e-11, 0.85, 6.478278331551995, 7.546431173856936, 0.0)
                         ,use_experiment=true);

# EEE EIS CV .. fitnes factos {1, 10} ... good fit for both ... only 2 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=100, bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true, 1.0, 1.0, 1.0, 24.383816208538917, 21.94866414517787, 20.68161992890677, 0.1380699169577035, -0.22933373067502555, -0.23743440380456354, 9.3e-11, 0.85, 2.7481343742171562, 2.930552174417333, 0.0)
                  ,use_experiment=true);                         
                         
### COOL ... without e_fac maybe the best EEE EIS CV fit                        
# EEE EIS CV ... fitness factors {1, 1} ... very good EIS fit ... 3 peaks     ######################## COOOL !!!!! ### but try pO2 = 80 !!!!                        
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=100, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true, 1.0, 1.0, 1.0, 21.96912506564684, 21.569512049276458, 20.93977987957153, 0.1241112223034025, -0.10660629726868907, -0.006580383278506099, 9.3e-11, 0.85, 5.848806563958082, 8.27404376722011, 0.0)
                  ,use_experiment=true);
                  
# EEE EIS CV ... fitness factors {0.2, 1} ... 2 peaks (+ 1 very small)      
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true, 1.0, 1.0, 1.0, 22.09624010786229, 21.66950610327097, 20.565105499265908, 0.1705724661650871, -0.1446619735003291, 0.034265126988040316, 9.3e-11, 0.85, 2.6704017916319613, 2.415004897727594, 0.0)
                  ,use_experiment=true);
                  
################## with e_fac ################# 

#!!!!!!! the following try fit with CV ... for pO2=80
# EEE EIS .... CV fitted automatically :) ... very good fit for .. with e_fac !!                  
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(true, 1.0, 1.0, 1.0, 21.9445, 21.3697, 20.8594, 0.101682, -0.178855, 0.115187, 9.3e-11, 0.85, 6.22922, 5.75616, 0.138391)
                         ,use_experiment=true);
                         
                         
                  

############################################################################
#################### LLL ###################################################
############################################################################


################# no_e_fac ##############
# 4 EIS ... not capturing the non-monotony ... even with extended search area ##### # # ## ## ## # # # 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC"     ], 
                                   prms_values=(0.0, 0.0, 0.0, 23.536793805295844, 26.44428245339161, 21.308492340939917, 0.5500221642340698, 0.029112401256050482, -0.020889037864824146, 9.3e-11, 0.85, 2.0593360762011694)
                                   ,use_experiment=true);                                        

# 4 EIS ELE ... example ... but very simmilar to"separate_vacancy", LLL                                 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                          prms_names["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "e_fac"], 
                                          prms_values=(0.0, 1.0, 0.0, 24.4631, 21.7136, 21.5347, -0.0113085, -0.0926479, 0.179689, 9.3e-11, 0.85, 3.70135, 0.0)
                                          ,use_experiment=true);
                                          
# 4 EIS EEL ... but pretty much the same as the above LEL                                          
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                          prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "e_fac"], 
                                          prms_values=(1.0, 1.0, 0.0, 24.1783, 21.7137, 21.5419, -0.0107187, -0.0862911, 0.185543, 9.3e-11, 0.85, 3.74939, 0.0)
                                          ,use_experiment=true);
                    
                    
################# with_e_fac ##############                    
                 
# EEL EIS .... CV hooked ... sensitive on pO2 ... 3 peaks... with e_fac
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(1.0, 1.0, 0.0, 21.9258, 21.2852, 21.7927, 0.0277078, -0.280821, 0.362658, 9.3e-11, 0.85, 7.54141, 1.65368, 0.296044)
                         ,use_experiment=true);
                         
# LEE EIS .... CV hooked
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(0.0, 1.0, 1.0, 24.253505704719995, 27.895308798638613, 20.689058014097917, 0.5435406241193925, 0.2124972805439068, -0.025572361570387683, 9.3e-11, 0.85, 2.1001141086263937, 0.8016033006905938, 0.004387336514558263)
                                ,use_experiment=true);

# ELE EIS .... CV hooked ... the worse sensitivity to pO2 ... and the worse shape of Nyquisty
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1.0, 0.0, 1.0, 22.004601154444064, 21.969771021412477, 20.64499552886603, -0.5553757496916619, -0.04999999999933978, 0.23448570623706932, 9.3e-11, 0.85, 5.925253670629901, 9.734714051365351, 0.7504653628307916)
                                ,use_experiment=true);                 
                 
#### LLL CV .... husty!!! ... trefa .... ale EIS je tragicke
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 24.2177, 21.6481, 22.6917, 0.698797, 0.252055, -0.122604, 9.3e-11, 0.85, 1.46538, 4.51906, 0.285749)
                                ,use_experiment=true);                    
##########
## COOL ##                                 
#### LLL EIS CV .... husty!!! ... obe nejsou spatni                                   
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(true, 0.0, 0.0, 0.0, 22.576839213753235, 22.011902293420093, 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, 0.3033398196955781)
                                              ,use_experiment=true);                    
                    
### LLL EIS CV .... BEZ e_fac   # --->> nic moc :)
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(true, 0.0, 0.0, 0.0, 24.069480130382026, 27.84830853588366, 22.089454852605634, 0.7999985765425789, -0.5103390865926676, -0.022909285494503932, 9.3e-11, 0.85, 0.9567488020824775, 2.851374822442069, 0.0)
                           ,use_experiment=true);                    
                    
                    
s

                                          
                                          
                                          
                                          
### HAND FIT ### :) LLL EIS ... quite good EIS... bad CV
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 
                                                     22.277, 22.312, 22.3, 
                                                     0.048, -0.271, 0.566, 
                                                     1.05e-11, 0.35, 24., 2.544, 0.303)
                                       ,use_experiment=true);
                                       
#                                                    #################################################
### corrected HAND FIT ... LLL EIS ... WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW WOWOWOWWOOW WO !!!!!! BRUTALISK !!!!! EIS cooool !!!!
#                                                    #################################################
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.325513444696398, 22.100276559559344, 22.382691084143286, 0.04243114428504254, -0.22279524777844634, 0.5430896356175956, 1.05e-11, 0.35, 34.397339021979015, 4.503994521917953, 0.338333333916094)
                                       ,use_experiment=true);
                                       
######### corrected HAND FIT for more pO2 ... LLL EIS 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true);
#### .... continuation from previos... to pO2 60
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.519825220050315, 22.101114333751948, 22.537917125481112, 0.046779285103424426, -0.036495661417036834, 0.6645252416469626, 1.06e-11, 0.35, 66.81244689470016, 2.899812097038787, 0.7357178655306754)
                                       ,use_experiment=true);
                           
##### ....... and ... then pO2 40, 60, 80,  again                           
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.49042860963425, 22.18758070122782, 22.535374543362902, 0.04626751671867485, -0.02797383209333109, 0.6910672027823139, 1.06e-11, 0.35, 51.10551528037439, 3.2209400174059213, 0.6775179360491119)
                                       ,use_experiment=true);
                    

#!#!!##!!#!#!#!#!#!#!#!#!#!#!#!#!!##!!#!#!#!#!#!#!##!#!!##!!#!#!#!#!#!#!#!#!#!#!#!#!!##!!#!#!#!#!#!#!#
### Temperature study !!#!#!#!#!!##!!#!#!#!#!#!#!#!#!#!#!#!#!!##!!#!#!#!#!#!#!#

#### only EIS fit
# 700 ... 2 EIS LLL ...     EIS very good even for pO2 20 !!! ... but CV bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(700, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.0248, 21.6486, 21.5078, -0.12619, 0.0824053, -0.129584, 1.02, 0.752653, 9.66799, 7.46362, 0.020354)
                                ,use_experiment=true); 
          
              # actaual good result
              new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(700, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.3641, 21.567, 21.3053, -0.304646, 0.405226, -0.266334, 1.02, 0.860279, 5.97275, 14.3188, 0.65)  
                                ,use_experiment=true);
                                
                                
                                
                          # NO _ EFAC
                          new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(700, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.0044, 21.6852, 21.4732, -0.00861625, 0.0414177, -0.0727471, 1.02, 0.850002, 8.59456, 7.56586, 0.0)
                                ,use_experiment=true);
                                
                                

                                                
                                
# 750 ... 2 EIS LLL ...    EIS great even for pO2 20 !!! ... but CV bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(750, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.3038, 22.0987, 22.0554, -0.209048, -0.0667363, 0.239006, 2.07, 0.38203, 27.384, 11.7518, 0.0148614)
                                ,use_experiment=true);


              # Actual good result
              new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(750, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false),  pyplot=1,  
              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
              prms_values=(1, 0.0, 0.0, 0.0, 22.4573, 21.8441, 21.9296, 0.238407, -0.395853, 0.348519, 2.07, 0.850001, 19.8925, 9.18612, 0.65)
                                ,use_experiment=true);
              
        ######    #  #####  # NO E_FAC !!!!              
                            new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(750, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                            prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                              prms_values=(1, 0.0, 0.0, 0.0, 22.2738, 22.0181, 21.9436, -0.0849979, -0.0602243, 0.168344, 2.07, 0.85, 16.6664, 9.38196, 0.0)
                            ,use_experiment=true);
          
              
              
                                
# 800 ... 2 EIS LLL ... 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true); 
    #alternative  
    new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 22.56783400557656, 22.22237966317968, 22.50693968907697, 0.05520337660564685, -0.033209183842112, 0.6770403513536281, 3.62, 0.35, 45.08895790653341, 3.269908550672186, 0.6182027814202923)
                                ,use_experiment=true);  
    
    # alternative #2
    new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 23.0467, 22.5346, 22.4665, -0.05, 0.0525283, 0.697041, 9.3e-11, 0.849408, 38.7172, 7.74472, 0.651215)
                                ,use_experiment=true); 
                                
    #alternative  #3 .. new one
    new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 23.671342687700914, 22.75009386895551, 22.36278147296789, -0.07877895777799904, 0.10251592543765028, 0.633939269456473, 3.62, 0.8503249682155306, 29.210203301013003, 12.194786069416887, 0.416244601465356)                                
                                ,use_experiment=true);
                                
            # Actual good result
            new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 23.075, 22.69, 22.4508, -0.21612, 0.0730355, 0.68861, 3.72, 0.850869, 37.1915, 16.7478, 0.65)                                
                                ,use_experiment=true);
            
#             # Actual good result #2 ->> ambigi-cosi--gious!! ... but behavior is the same... CV and EIS for [20, 40, 60, 80] ... strange !!!!
#             new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
#                                 prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
#                                  prms_values=(1, 0.0, 0.0, 0.0, 23.8216, 22.9161, 22.488, -0.0877746, 0.306968, 0.704817, 3.72, 0.850195, 41.096, 49.1646, 0.65)
#                                 ,use_experiment=true);                              
                 
  ######    #  #####  # NO E_FAC !!!!
                      new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.7906, 22.3002, 21.8264, -0.288601, 0.134113, 0.085357, 3.72, 0.850009, 7.11877, 28.32, 0.0)
                                ,use_experiment=true);                 
                 
                 
# 850 ... 2 EIS LLL ... 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 25.474103553624037, 22.925019737001634, 21.85804785753094, -0.10831702948372236, 0.38984727556707693, -0.23477140773747657, 5.85, 0.65, 3.4114628080433516, 20.67717563382255, 0.0)
                                       ,use_experiment=true);
                                       
  ######    # Actual good result #2 ->> ambigi-cosi--gious!! ... but behavior is the same... CV and EIS for [20, 40, 60, 80] ... strange !!!!
            new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 25.0252, 22.1444, 21.6018, -0.0807316, 0.400022, -0.224248, 5.85, 0.89, 2.19752, 2.81139, 0.65)
                                ,use_experiment=true); 
             
            
  ######    #  #####  # NO E_FAC !!!!
                      new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 24.4626, 22.6312, 21.9401, 0.338431, 0.0088969, 0.169719, 5.85, 0.850069, 3.95899, 6.81854, 0.0)
                                ,use_experiment=true);                                                                                                                                
                                   
  ######    #  #####  # NO E_FAC !!!! pO2 = 40, 60
                      new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 25.2697, 22.6163, 21.8379, -0.17644, 0.213803, -0.209637, 5.85, 0.895986, 3.26162, 5.97785, 0.0)
                                ,use_experiment=true);                                
                              
                              
                                   
                                   
           
#### EIS-CV fits


## 700 ... 2 EIS-CV LLL  ... ff1-100 --> ok .. :) we can take it
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=700, pO2=[20, 40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(1, 0.0, 0.0, 0.0, 
                                                            21.9222, 21.2101, 21.809,
                                                            0.0558328, 0.0480154, -0.0438488, 
                                                            1.02, 0.765282, 37.9658, 5.55807, 0.270473)
                                              ,use_experiment=true);
           
## 750 ... 3 EIS-CV LLL ... ff ???
# first attempt with ff1-100 --> EIS too bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(750, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(1, 0.0, 0.0, 0.0, 22.2922, 21.4767, 26.0711, 0.0326867, -0.0220428, 0.0672773, 2.07, 0.83353, 34.7543, 7.64352, 0.251502)
                                              ,use_experiment=true);
 
# second attempt with ff5-100 --> much better, 3 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=750, pO2=[40, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(1, 0.0, 0.0, 0.0, 
                                                           22.178, 21.508, 22.1441, 
                                                           0.00928035, 0.210568, -0.104948, 
                                                           2.07, 0.475152, 41.6055, 9.92313, 0.339722)
                                                     ,use_experiment=true);

# third attempt with ff10-100 -->  ok .. :) we can take it
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=750, pO2=[40, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(1, 0.0, 0.0, 0.0, 
                                                          22.1636, 21.4919, 22.0697,  
                                                          -0.00599031, 0.18108, -0.05,  
                                                          2.07, 0.532246, 51.0891, 9.42743, 0.362266)
                                                     ,use_experiment=true);                                                     
                                                     

## 800 ... 2 EIS_CV LLL ... ff ???
# some previous attempt ... unknown ff ... only 2 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(true, 0.0, 0.0, 0.0, 
                                               22.576839213753235, 22.011902293420093, 27.69994335225577, 
                                               0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 
                                               9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, 0.3033398196955781)
                                              ,use_experiment=true);

# some result .. .. 
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(0.0, 0.0, 0.0, 22.485999537854475, 22.17094767696965, 25.601561532601615, 0.027450455830081023, -0.35838934449884186, 0.7546179586337105, 9.3e-11, 0.85, 55.73516673393179, 9.18503138670932, 0.33473609546853517)
                                              ,use_experiment=true);
                    
# sedonc???
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(1, 0.0, 0.0, 0.0, 22.4555, 21.6669, 22.7248, -0.0285608, 0.195182, 0.102502, 3.58, 0.0502756, 28.2008, 17.8957, 0.4)
                                              ,use_experiment=true);





## 850 ... 2 EIS_CV LLL ... ff ???
           
           
           
           
           
           
###########################################################################
# Summary of temperature study:
################  only EIS ################################################
# pO2 = 40, 60
#                          "A.r",   "R.r",   "O.r",    "A.DG",    "R.DG",     "O.DG",     "conductivity", "nu",      "OC",    "ms_par",  "e_fac"]
#                         
# 700: (1, 0.0, 0.0, 0.0, 22.3641,  21.567,  21.3053, -0.304646,   0.405226,  -0.266334,      1.02,       0.86,      5.97275, 14.3188,    0.65)
# 750: (1, 0.0, 0.0, 0.0, 22.4573,  21.8441, 21.9296,  0.238407,  -0.395853,   0.348519,      2.07,       0.85,     19.8925,   9.18612,   0.65)
# 800: (1, 0.0, 0.0, 0.0, 23.075,   22.69,   22.4508, -0.21612,    0.0730355,  0.68861,       3.72,       0.85,     37.1915,  16.7478,    0.65)
# 850: (1, 0.0, 0.0, 0.0, 25.0252,  22.1444, 21.6018, -0.0807316,  0.400022,  -0.224248,      5.85,       0.89,      2.19752,  2.81139,   0.65)
################
# pO2 = 20, 40, 60
# zeta = 0                 
#                         
# 700: (1, 0.0, 0.0, 0.0, 22.0044,  21.6852, 21.4732, -0.00861625, 0.0414177, -0.0727471,     1.02,       0.850002,  8.59456,  7.56586,   0.0)
# 750: (1, 0.0, 0.0, 0.0, 22.2738,  22.0181, 21.9436, -0.0849979, -0.0602243,  0.168344,      2.07,       0.85,     16.6664,   9.38196,   0.0)
# 800: (1, 0.0, 0.0, 0.0, 22.7906,  22.3002, 21.8264, -0.288601,   0.134113,   0.085357,      3.72,       0.850009,  7.11877, 28.32,      0.0)
# 850: (1, 0.0, 0.0, 0.0, 24.4626,  22.6312, 21.9401,  0.338431,   0.0088969,  0.169719,      5.85,       0.850069,  3.95899,  6.81854,   0.0)
################
# pO2 = 40, 60, zeta=0
#
# 700: (1, 0.0, 0.0, 0.0, 22.0245, 21.7055, 21.5085,  -0.0852986,  0.0959023, -0.16008,       1.02,       0.850003,  8.96497,  7.57717,   0.0)
# 750: (1, 0.0, 0.0, 0.0, 22.2717, 22.1757, 22.0818,   0.0421601, -0.222154,   0.312441,      2.07,       0.850457, 23.5146,   7.42426,   0.0)
# 800: (1, 0.0, 0.0, 0.0, 23.3087, 23.1146, 22.4228,  -0.0944925, -0.0365465,  0.59062,       3.72,       0.871189, 28.9727,  22.5948,    0.0)
# 850: (1, 0.0, 0.0, 0.0, 25.2697, 22.6163, 21.8379,  -0.17644,    0.213803,  -0.209637,      5.85,       0.895986,  3.26162,  5.97785,   0.0)
#
################  EIS-CV   ################################################
#                          "A.r",   "R.r",   "O.r",    "A.DG",    "R.DG",     "O.DG",     "conductivity", "nu",   "OC",    "ms_par",  "e_fac"]
#                         
# 700: (1, 0.0, 0.0, 0.0, 21.9222, 21.2101, 21.809,    0.0558328,  0.0480154, -0.0438488,     1.02,       0.765282, 37.9658,   5.55807,   0.270473)
# 750: (1, 0.0, 0.0, 0.0, 22.1636, 21.4919, 22.0697,  -0.00599031, 0.18108,   -0.05,          2.07,       0.532246, 51.0891,   9.42743,   0.362266)
# 800:  ???
# 850:  ???
###########################################################################

#### Temperature initial condition inspired from the above
#["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
#                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
#                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
#                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
#                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"]#      
# Not BAD !!! :) :) :) for both EIS and CV 
        (1, 0.32, 0, 0, 0,
        0.2414, 21.9222,       -0.0618,     0.0558328, 
        0.2818, 21.2101,       0.133,     0.0480154,
        0.2607, 21.809,        -0.00615,    -0.0438488,
        -0.133,  0.765282,     13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807)   
#       
#       
#       
#       (0, 0.0, 0.0, 0.0, 0.24139999999999873, 0.2818000000000005, 0.26069999999999993, -0.06182311, 0.13306459999999998, -0.0061512000000000025, 1.0499999999999998, -0.23303600000000002, 13.1233, 3.8693599999999995, 0.09179299999999996)
############################  


### #
#
###### 

 # 3  # 3 # #
  #                                # IDEA: Take multiple "the best" fits for every condition and make many "probability curves"
  
   #  # ##K # #L KJ# L#K JL 
           
           
           
           
# fit pres bety a S ... nic moc stejne! :( 
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=700, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1, 
                                prms_names=["separate_vacancy",
              "A.beta", "R.beta", "O.beta",
              "A.S", "R.S", "O.S",
              "A.exp", "R.exp", "O.exp", 
              "A.r", "R.r", "O.r",              
              "A.DG", "R.DG", "O.DG",     
              "conductivity", "nu",      "OC", "ms_par", "e_fac"  ],
                                 prms_values=(1, 0.5468986658625657, 0.47998951044726024, 0.2779139613165412, -0.27769563747657716, -0.09593075662029801, -0.22320858430182988, 0.0, 0.0, 0.0, 21.9464, 21.5159, 21.4928, 0.0032943, 0.0277383, -0.0838499, 1.02, 0.811907, 10.472, 6.00553, 0.0718344)
                                ,use_experiment=true);             
           
           
           
           
           
           
           
           
           
           

           

           
           
           
#######################################################################################################################################
################################################# shared ##############################################################################
#######################################################################################################################################

############################################################################
#################### EEE ###################################################
############################################################################

                                          
##### MAYBEE...... another hand-fit 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 
                                                            21.9, 21.3, 21.0, 
                                                            0.15, 0.0, -0.40967, 
                                                            1.06e-11, 0.35, 1.0, 20., 0.0)
                                                ,use_experiment=true);
############
### COOL ###
### EEE 1 EIS ... sha ... hustYYY!! ...                                                 
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 21.78592822930247, 21.347238180487707, 20.845011931720087, 0.11415904394451544, -0.04710867464932109, -0.4221231346344268, 1.06e-11, 0.35, 1.252221527673722, 37.34152581968751, 0.0)
                                                ,use_experiment=true);                                                
                                                
                           
### EEE 1 EIS-CV ... sha ... fitness = {1 ,1}   --->> great EIS, medium CV
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 21.77097654529923, 21.339094620305584, 20.83825101342013, 0.11612661311085312, -0.01621306293886935, -0.44405718790607684, 1.06e-11, 0.35, 0.6404908464083954, 37.19354903267066, 0.010059997045145626) 
                                                ,use_experiment=true);
                           
### EEE 1 EIS-CV ... sha ... fitness = {0.5 ,1} --->> very good EIS, not good CV                          
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 21.768338176697593, 21.339959521197464, 20.77053224186594, 0.14211354657568556, -0.06522415408171552, -0.4099634528362534, 1.06e-11, 0.35, 0.9838424541352313, 43.614438009170286, 0.04296470330868719)
                                                ,use_experiment=true);                           
                           
### EEE 1 EIS-CV ... sha ... fitness = {0.1 ,1} --->> EIS 2 peaks, good CV                         
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 27.204738430905046, 21.042465727359954, 20.903492276526453, -0.7650977954784123, -0.16483434129293068, -0.5474897673095572, 1.06e-11, 0.35, 7.553853360574966, 20.131456245497063, 0.12330883448476397)
                                                ,use_experiment=true);                            
                           
### EEE EIS CV.... shared vac ... not bad, not the best ... 2 peaks
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                          prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC",  "ms_par","e_fac"     ], 
                                          prms_values=(false, 1.0, 1.0, 1.0, 21.6797, 21.0495, 20.7467, 0.0212226, -0.423638, 0.798779, 9.3e-11, 0.85, 6.75588, 9.95628, 0.305956)
                                          ,use_experiment=true);                                        
                          
### EEE 1 CV ... sha ...   --->> EIS 1 peak, better CV                         
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 27.280316820837967, 20.744018794803896, 25.703864047655617, -0.10168151776708975, -0.10239784764049484, -0.7893449634234231, 1.06e-11, 0.35, 0.0029974884900697343, 8.174364689405985, 0.11692539726030851)                           
                                                ,use_experiment=true);                                           
                                          

############################################################################
#################### LLL ###################################################
############################################################################


#### LLL 1 EIS .... shared vac   -> CV bad , EIS 2 peaks

#### ...... LLL shared ... 1 EIS CV ... fitness = {10, 10}
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(0, 0.0, 0.0, 0.0, 22.6529, 21.2666, 22.653, 0.201058, 0.0140769, -0.223763, 1.06e-11, 0.35, 1.33437, 25.0502, 0.644818)
                                              ,use_experiment=true);
                                              
#### LLL shared ... 1 EIS CV .. fitness = {1, 10}
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(false, 0.0, 0.0, 0.0, 23.660999804374093, 21.16644727635489, 23.246774488288892, 0.5590363181503285, 0.3057035105534315, -0.46167785857357707, 1.06e-11, 0.35, 57.20049929904167, 10.949218885809483, 0.40238189628810334)
                                                     ,use_experiment=true);                                              

#### LLL shared ... 1 EIS CV .. fitness = {1, 5}                                                     
ysz_fitting.simple_run(simulations=["EIS", "CV"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                       prms_values=(0, 0.0, 0.0, 0.0, 22.3718, 22.159, 26.9917, 0.00119476, -0.217313, 0.74843, 1.06e-11, 0.35, 12.344, 11.4685, 0.383988)
                                                       
                                                       ,use_experiment=true);                                                     

#### ...... LLL shared !!! 2 EIS CV ...not bad :)   >>>>>>> almost THE SAME result as for the separate vac
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(false, 0.0, 0.0, 0.0, 22.5615, 22.113, 27.0515, 0.0456777, -0.33161, 0.684794, 9.3e-11, 0.85, 75.0419, 9.5712, 0.311755)
                                              ,use_experiment=true);



                                                                         

                           
                           
                           
                           
                           
                           
                           
                           
        
                           
                           
                            

                           

                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
###########################################################                           
#####################  PAR_STUDY  #########################                           
# separate VS shared vacancy ... 3 vs 2 peaks
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(true, 1.0, 1.0, 1.0, 
                                       21.6106, 21.898, 21.4851, 
                                       0.497806, collect(-0.3 : -0.02 : -0.3), 0.229967, 
                                       9.3e-11, 0.85, collect(1.0 : 1 : 14), 3.5, 0.0)
                           ,use_experiment=false);

                           
# e_fac descussion ... without inductance
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(0.0, true, 0.0, 0.0, 0.0, 
                                       23.6106, 22.598, 23.4851, 
                                       -0.497806, collect(-0.3 : -0.02 : -0.3), 0.229967, 
                                       9.3e-11, 0.85, 1.0, 1.0, collect(0.0 : 0.1 : 0.5))
                           ,use_experiment=false);
                           
                           
                           
                           
                           
                           
# OC study ... 3 clear peaks  LLL   !!!!                       
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(0.0, [true, false], 0.0, 0.0, 0.0, 
                                       22.2106, 22.098, 21.5851, 
                                       -0.1497806, 0.0, 0.029967, 
                                       9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.1 : 0.0))
                           ,use_experiment=false); 
# separate_vacancy study ... 3 clear peaks ... LLL !!!
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, [true, false], 0.0, 0.0, 0.0, 
                                              22.2106, 21.798, 21.1851, 
                                              -0.1497806, -0.0, -0.19967, 
                                              9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.1 : 0.0))
                                  ,use_experiment=false);
                           
# OC study ... 3 clear peaks ... EEE !!!!                          
                           ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                              22.1106, 22.098, 21.5851, 
                                              0.1497806, -0.0, 0.029967, 
                                              9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                  ,use_experiment=false);
                           
                           
########################################################################################################                           
#  .... .. .. 4 peaks clear!... EEE !!! bulk peak!
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                       prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                       prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                                   21.2, 21., 21.0, 
                                                                   0., 0., 0., 
                                                                   [0.05]*1.0e-11 , [0.9], 30., 10., collect(0.0 : 0.1 : 0.0))
                                                       ,use_experiment=false);
########################################################################################################
                           
                           
                           
# e_fac study .... 3 clear peaks ... down loping !!!! :D
# the same looping is appearing using EXP kinetics. 
#### when A or O is fast, middle peak starts to be the boundary one.... and makes a "inductance tail" 
# --- dva nebo tri obrazky
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, true, 0.0, 0.0, 0.0, 
                                              22.2106, 22.098, 21.5851, 
                                              -0.1497806, 0.0, 0.029967, 
                                              9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.2 : 0.5))
                                  ,use_experiment=false);


                                  
### DD vs nu  --->> DD changes only R_ohm ... nu changes DRT
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                            22.0, 22., 21.5, 
                                                            0.15, 0.0, -0.629967, 
                                                            [0.1, 1, 10, 100].*1.0e-11, collect(0.1 : 666.1 : 0.8), 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false); 
                                  
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
###########################################################################################
##########################        USED IN ARTICLE      ####################################
###########################################################################################
###########################################################################################
###########################################################################################

# template for saving files
TC=[700, 750, 800, 850]; pO2=[20, 40, 60]; data_set="OLD_MONO_100";
       ysz_fitting.simple_run([
          ysz_fitting.CV_simulation(TC, pO2, data_set=data_set, sample=50, physical_model_name="ysz_model_GAS_LoMA_generic"  )...,
          ysz_fitting.EIS_simulation(TC, pO2, data_set=data_set, f_range=(1, 30000, 1.2), physical_model_name="ysz_model_GAS_LoMA_generic"  )...],
                                                      save_files=true, save_dir="experiments",  pyplot=0,  
                                                      prms_names=Nothing,
                                                      prms_values=[]
                                                      ,use_experiment=true, line_color_idx=1);


# par_stud_e_fac_LF
new_GAS_LoMA_shared = ysz_fitting.simple_run(
                                                       ysz_fitting.EIS_simulation(800, [80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                       #ysz_fitting.CV_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                       
                                                           pyplot=1,  
                                                            prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                             prms_values=(0.0, 0.0, 0.0, 22.5, collect(22.9 : 0.1 : 22.9), 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, collect(0.0 : 0.05 : 0.2))
                                                            ,use_experiment=false);
                                                            
# par_stud_e_fac_HF
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                          prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"], 
                                          prms_values=(1.0, 1.0, 1.0, 21.6106, 27.898, 20.4851, -0.497806, -0.234529, 0.229967, 9.3e-11, 0.85, 3.04504, 0.05, [0.0, 0.1, 0.2, 0.3, 0.4])
                                          ,use_experiment=false);

                                          
# par_study_sha_sep_2_peaks >>> LaTeX tab. parameters: id 0.1
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                         prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                         prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                     22.0, 22., 21.5, 
                                                     0., 0.0, 0., 
                                                     9.3e-11, 0.5, 300., 100., collect(0.0 : 0.1 : 0.0))
                                         ,use_experiment=false);

# par_study_sha_sep_3_peaks >>> LaTeX tab. parameters: id 0.2
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                         prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                         prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                     22.0, 22., 21.5, 
                                                     0.15, 0.0, -0.629967, 
                                                     9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                         ,use_experiment=false);

# VERY interesting !!!
# separare_vacancy interchanged the order of reactions <-> peaks bijections >>> LaTeX tab.: id 0.2
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                            22.0, [22., 29.], 21.5, 
                                                            0.15, 0.0, -0.629967, 
                                                            9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false);

# par_study_DGR >>> LaTeX tab. parameters: id 0.1                                         
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                            22.0, 22., 21.5, 
                                                            0., collect(-0.3 : 0.1 : 0.1), 0., 
                                                            9.3e-11, 0.5, 300., 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false);
                                                
# HOT candidate !!! corrected_hand_guess -> expanded to 3 pO2 ..... CV rubbish
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true);

                                       
#### SHA vs SEP comparison ##########
# shared
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(0, 0.0, 0.0, 0.0, 22.6529, 21.2666, 22.653, 0.201058, 0.0140769, -0.223763, 1.06e-11, 0.35, 1.33437, 25.0502, 0.644818)
                                                     ,use_experiment=true);
        
# separater        
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(1, 0.0, 0.0, 0.0, 22.325513444696398, 22.100276559559344, 22.382691084143286, 0.04243114428504254, -0.22279524777844634, 0.5430896356175956, 1.05e-11, 0.35, 34.397339021979015, 4.503994521917953, 0.338333333916094)
                                              ,use_experiment=true);

#### EXP vs LoMA comparison ########
# EEE 1 EIS-CV ... sha ... fitness = {0.1 ,1} --->> EIS 2 peaks, good CV                         
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 27.204738430905046, 21.042465727359954, 20.903492276526453, -0.7650977954784123, -0.16483434129293068, -0.5474897673095572, 1.06e-11, 0.35, 7.553853360574966, 20.131456245497063, 0.12330883448476397)
                                                ,use_experiment=true);
                                                
# LLL 1 EIS-CV ...fail of CV                       
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp ", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "conductivity", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(true, 0.0, 0.0, 0.0, 22.790631052263787, 22.30016049493873, 21.826412050023904, -0.2886011047255152, 0.13411306312611201, 0.08535699920265348, 3.72, 0.8500086485807636, 7.118766097564314, 28.32000294249703, 0.0)
                                                ,use_experiment=true);
                                                
# good LoMA for both ... 800 , [20, 60]                                                
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation([ 800], [20, 60], 0.0, f_range=(1. : 30000 : 1.2), physical_model_name="ysz_model_GAS_LoMA_Temperature", data_set="OLD_MONO_100", plot_legend=false), 
                      pyplot=1,   
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.4, 0.0, 0.0, 0.0, 
                      0.239, 21.956, -0.004, 0.073, 0.278, 21.096, -0.056, 0.08, 0.731, 21.823, 0.151, -0.043, -0.254,
                      
                      0.76,
                      
                      9.985, 51.946, 0.19, 0.797)
                      ,use_experiment=true);

                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
###########################################################
###########################################################
#### SHA vs SEP comparison ##########    CORRECT MODEL #### 
# shared
SHA_vs_SEP = ysz_fitting.simple_run(simulations=["EIS"], TC=750, pO2=60, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                        physical_model_name="ysz_model_GAS_LoMA_generic",          
                        prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r",      "A.DG",
                                        "R.r",      "R.DG",
                                        "O.r",      "O.DG",
                                        "nu",     "CO",     "COmm"],
                        prms_values=(0, 0.0, 0, 0, 0, 22.986833011976863, 0.09356569189369754, 21.955873372516326, 0.07492013587080583, [22.67, 29], -0.6154237534608018, 0.3414399051362852, 666, 17.26576136395294)
                        ,use_experiment=true);                     
                      
# separate
SHA_vs_SEP = ysz_fitting.simple_run(simulations=["EIS"], TC=750, pO2=60, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                        physical_model_name="ysz_model_GAS_LoMA_generic",          
                        prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                        "A.r",      "A.DG",
                                        "R.r",      "R.DG",
                                        "O.r",      "O.DG",
                                        "nu",     "CO",     "COmm"],
                        prms_values=(1, 0.0, 0, 0, 0, 22.239645778031775, -0.11475291933274577, 22.01350710083606, 0.08185420233832269, [21.95, 29], -0.21962624735381847, 0.3249340342508185, 21.975075442182572, 1.7008124346139954)
                        ,use_experiment=true);                     

#### EXP vs LoMA comparison ########
# EIS-CV ff1:1 ...          >>> OK                   
EXP_vs_LOMA = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=800, pO2=[20, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                        physical_model_name="ysz_model_GAS_LoMA_generic",          
                        prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                    "A.r",      "A.DG",   "A.S",
                                    "R.r",      "R.DG",   "R.S",
                                    "O.r",      "O.DG",   "O.S",
                                    "nu",     "CO",     "COmm"],
                        prms_values=(1, 0.17, 1, 1, 1, 21.9, 0.161, -0.416, 21.23, -0.126, 0.58, 20.977, 0.071, -0.005, 0.42, 14.052, 1.203)
                        ,use_experiment=true);  
                                                
# EIS-CV ff1:1  ... 3 peaks & fail of CV    >>> OK             
EXP_vs_LOMA = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=800, pO2=[20, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                        physical_model_name="ysz_model_GAS_LoMA_generic",          
                        prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                    "A.r",      "A.DG",   "A.S",
                                    "R.r",      "R.DG",   "R.S",
                                    "O.r",      "O.DG",   "O.S",
                                    "nu",     "CO",     "COmm"],
                        prms_values=(1, 0.0, 0, 0, 0, 22.052, -0.043, -0.511, 23.164, 0.114, 0.363, 22.001, 0.014, 0.019, 0.3, 12.655, 3.597)
                        ,use_experiment=true);
                        
# EIS-CV ff1:1  ... 2 peaks & CV ok         >>> OK             
EXP_vs_LOMA = ysz_fitting.simple_run(simulations=["EIS", "CV(f)"], TC=800, pO2=[20, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                        physical_model_name="ysz_model_GAS_LoMA_generic",          
                        prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                                    "A.r",      "A.DG",   "A.S",
                                    "R.r",      "R.DG",   "R.S",
                                    "O.r",      "O.DG",   "O.S",
                                    "nu",     "CO",     "COmm"],
                        prms_values=(1, 0.8, 0, 0, 0, 22.227, 0.005, -0.387, 22.123, 0.25, 0.327, 24.989, -0.013, 0.016, 0.687, 20.188, 1.38)
                        ,use_experiment=true);                        

                      
##### por-0 ##########

                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
# comparison of zeta VS LSM

prms_names=["separate_vacancy",         "A.r",     "A.DG",     "A.S",
                                        "R.r",     "R.DG",     "R.S",
                                        "O.r",     "O.DG",     "O.S",
                                        "nu",    "CO",     "COmm"],
                                      
prms_values=(1, 21.517,  -0.094, -0.33, 
                21.732,  0.202, 0.33, 
                 21.215,  -0.156, 0.0, 
                 0.75,       4.116,  0.38)   
                      
                      
                      
                      
                      
                      
                      
