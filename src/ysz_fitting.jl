module ysz_fitting
################################################################################################
####### TODO ###################################################################################
# [x] better ramp ... starting directly from steadystate :)
# [x] general global search using projection to each variable
# [x] compute EIS exacly on checknodes and therefore remove plenty of "EIS_apply_checknodes"
# [x] put appropriate and finished stuff into "CV_fitting_supporting_stuff"
# [x] postarat se o modularitu ysz_model_COSI, at se nemusi vytvaret znovu "experiments_COSI"
# [x] vymyslet lepe zadavani parametru pro ruzne modely s ruznymi parametry 
# [x] opravdu promyslet to objektove programovatni ... CV_sim, EIS_sim ... prms_lists ...
# ---[x] pridat treba dalsi experiment s dalsimi daty, vuci kterym se da srovnavat (kapacitance)
# ---[x] udelat tridu par_study
# [x] snehurkove fitovani by slo zrychlit, kdyz bych skriptoval EQ parametry a job by menil jen kineticke? ... chrm ...
# ---[x] nejak si preposilat steadystate?! 
# ------[x] zacina se ze steady_statu, ktery se pocita rychle
# [N] simple_run by mohl vracet par_study
# ---[N] par study by pak mohla mit funcki "uloz me"
# [N] prms_names a prms_values dat do jedne tridy
# ---[N] vlastne jakoby nevim, jestli je to dobry napad?
# 
# [ ] do experiments.jl pridat obecne zaznamenavani promennych od final_plot
# [ ] f_range pouzivat jako Array frekvenci, nikoliv jako trojici cisel
# [ ] konzistentne pridat "f_interval" do vsech EIS
# [ ] zaradit jednoduchy vztah teplotni zavislosti na odporu
# [x] dodelat v CV simulaci vse potrebne, co je v EIS. I importovani novych souboru.
# [ ] f_interval pridat do SIM_fitting
# [ ] fitness funkce s maximovou metrikou
# [ ] zobrazovani plot_EEC_data_general(...) pro "R3 + R4"
# [!!!] finally SEPARATE par_study and SIM_fitting meta function
# [ ] lepe vymyslet zadavani v externim souboru SIM_fitting
# [ ] plot slurm results -> spravne stridat barvy & styly cary
#
#
#
#
#
# ##### snehurka stuff ############################
# [x] snehurka, velke objemy dat, maska a metafile
# ---[x] save_dir preposilany parametrem a odlisit scripted_dir
# ---[x] sneh - zadavani i jinych parametru nez prms
# ---[x] sneh - pouzivat jen jeden soubor pro spouteni sbatch ... v hlavicce #!(..)/julia
# ---[N] mozna pouzit precompile
# ---[!!!!!] nechat meta_run_par_study() vytvorit skriptovaci soubor (dle poctu parametru, module, pocet experimentu)
# 
#
#
#
####### TODO - less important ###################################################################################
#
# ##### I-V stacionarni krivky
# [ ] zatim to modeluji jako CV s nizkym voltratem
# [ ] make I-V simulation
#
# ##### CAP simulation ############################
# [ ] add data_set and other features to CAP 
#
# ##### par-study #################################
# [x] PAR_STUDY vyhodnoceni ... soubory do slozky dane studie .. automatizovane
# ---[ ] pri zobrazovani dat udelat clustery dle chyby/prms a pak zobrazit (prms ci Nyquist) 1 reprezentanta z kazde
# ---[ ] automatizovat ukladani souboru vysledku par_study.vyhodnoceni() 
# [ ] promyslet velikost Floatu pri importu velkeho mnozstvi dat
# [x] funkci pro par_study, aby umela zmenit sve info (redukovat pO2_list a tak) -> funkce par_study_filter()
# [ ] par_study_plot_the_best() by mela vyhodit jeden figure s dvema subploty... asi ... 
# [ ] snehurka by rada ./zabal.sh a pocitac zase ./rozbal_par_study.sh
# [ ] automaticke vybirani vhodne scripted_tuple
# 
# ##### interpolacni fitovani #####################
# [x] vizualizace trendu chyby mezi interpolanty
# ---[x] zobrazovat soucty odchylek
# ------[x] obrazky po normalizaci vypadaji ruznorode
# ---[x] zkusit vykoukat trend z nenormalizovanych dat
# [?] aplikovat masku pro fitting pro scan_2D_recursive
#
#
# #### Fitting process ######
# [x] prozkoumat experimentalni data (a udelat prislusne procedury)
# [x] ysz_model_GAS_LoMA_Temperature
# ---[x] find good fits for each TC -> 3 peaks, reasonable CV -> at least 700, 750
# ---[ ] perform TEMP fitting with this initial guess
#                   
###################################################################################################



using Printf
using PyPlot
#using Plots
using DataFrames
using LeastSquaresOptim
#using Optim
using BlackBoxOptim
using LinearAlgebra

import Base.string



#include("../examples/ysz_experiments.jl")


include("../src/general_supporting_stuff.jl")
include("../src/import_experimental_data.jl")
include("../src/export_simulated_data.jl")

include("../src/DRT.jl")
include("../src/EIS_preprocessing.jl")

include("../src/simulations/general_simulation.jl")
include("../src/simulations/CV_simulation.jl")
include("../src/simulations/EIS_simulation.jl")
include("../src/simulations/CAP_simulation.jl")

include("../src/EIS_preprocessing_tool.jl")

include("../src/par_study.jl")
#include("../src/interpolation_fitting.jl")
include("../src/EEC_module.jl")


function construct_SIM_list(;TC=Nothing, pO2=Nothing, bias=0.0, data_set=Nothing, simulations=Array{String}(undef, 0),
                    fitness_factors=Nothing, physical_model_name=Nothing)
  SIM_list = []
  try
    if fitness_factors == Nothing
      aux_array = zeros(length(simulations))
      aux_array .= 1
      fitness_factors = aux_array
    end
    SIM_list = get_SIM_list_rectangle(TC, pO2, bias, data_set, simulations, fitness_factors, physical_model_name)
  catch e
    println(e)
    println("ERROR: please define TC, pO2, bias, data_set, simulations ... OR ... define SIM_fitting")
    return throw(Exception)
  end
  return SIM_list
end

function simple_run(SIM_list=Nothing; TC=Nothing, pO2=Nothing, bias=0.0, data_set=Nothing, simulations=Array{String}(undef, 0),
                    fitness_factors=Nothing, physical_model_name=Nothing, 
                        pyplot=0, use_experiment=true, prms_names=[], prms_values=[], 
                        test=false, save_files=false, save_dir="default",
                        line_color_idx = 1, evaluate_temp=false
                        )
  # here starts the true body
  if SIM_list==Nothing
    SIM_list = construct_SIM_list(TC=TC, pO2=pO2, bias=bias, data_set=data_set, simulations=simulations, 
                                    fitness_factors=fitness_factors, physical_model_name=physical_model_name)
  end
  
  save_path="../data/simple_run/"*save_dir
  
  if evaluate_temp
    plot_temp_parameters(prms_values=prms_values, prms_names=prms_names, label=save_dir, 
                          save_file=save_files, save_dir=save_path*"/", file_name="temp_prms",
                          line_color_idx=line_color_idx)  
  end
  
  if test
    test_result = 0
  end
  res = DataFrame()
  for SIM in SIM_list
    
    
    function recursive_simple_run_call(output_prms, plot_names, plot_values, active_idx)
      if active_idx > size(prms_names,1)

        ###################################################################
        # for each combination of output_prms from prms_lists perform this #######
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, output_prms)
      
        
        if check_equal_size(prms_names, prms_values)  
          SIM_raw = typical_run_simulation(SIM, prms_names_in, prms_values_in, pyplot)        
        end    
        if size(plot_names,1) < 1
          plot_prms_string = ""
        else
          plot_prms_string = " $(string(plot_names)) = $(plot_values)"          
        end
        SIM_sim = apply_checknodes(SIM, SIM_raw, SIM.checknodes)
        
        ## TODO !!!
        res = SIM_sim
        ###
        
        if pyplot > 0
            typical_plot_sim(SIM, SIM_sim, plot_prms_string)
        end
        if save_files
          save_file_prms(SIM, SIM_sim, save_path, plot_values, plot_names, [], mode="sim")
        end
        if use_experiment
          if test
            test_result+=fitnessFunction(SIM, SIM_sim, SIM_exp)
          else
            fitness_error_report(SIM, plot_prms_string, SIM_exp, SIM_sim)
          end
        end
        return
        ###################################################################
        ###################################################################
      
      end
      if size(prms_values[active_idx],1)>1
        for i in prms_values[active_idx]          
          recursive_simple_run_call(
            push!(deepcopy(output_prms),i),
            push!(deepcopy(plot_names),prms_names[active_idx]),
            push!(deepcopy(plot_values),i),
            active_idx + 1)
        end
      else
        recursive_simple_run_call(
          push!(deepcopy(output_prms),prms_values[active_idx][1]),
          plot_names,
          plot_values,
          active_idx + 1)
      end
    end

    # here the true body continues
    if use_experiment
      SIM_exp = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
      if (pyplot > 0)
        typical_plot_exp(SIM, SIM_exp)
      end
      if save_files
        save_file_prms(SIM, SIM_exp, save_path, [], [], [], mode="exp")
      end
    end
    if prms_names!=Nothing
      recursive_simple_run_call([], Array{String}(undef,(0)), Array{Any}(undef,(0)), 1)
    end
  end
  
  if test
    return test_result
  else
    return res
  end
end



###########################################################
###########################################################
#### Fitness_function_interpolation #######################
###########################################################
###########################################################




###########################################################
###########################################################
####     D R T     ########################################
###########################################################
###########################################################



function test_DRT(SIM_list=Nothing;lambda=0.0, mode="EEC", TC=800, pO2=80, bias=0.0, R_ohm=1, R1=1, C1=0.001, R2=1, C2=0.0001, alpha=1, prms_names=[], prms_values=[], backward_check=true, draw_semicircles=false, plot_option="Nyq DRT Bode RC", f_range=EIS_get_shared_f_range(), data_set="MONO_110", physical_model_name=Nothing,
tau_min_fac=10, tau_max_fac=10, tau_range_fac=2,
peak_merge_tol=0.0, fig_num=EIS_standard_figure_num, plot_bool=true,
CAP_comparison=false, CAP_bottleneck_prm="rR", CAP_plot_CAP_CV=true, CAP_plot_num = 101)
  

  if data_set=="POLY_OCV_test"
    bias=0.0
    data_set_list = ["POLY_$idx" for idx in [1,2,3]]
  else  
    data_set_list = [data_set]
  end
  
  if CAP_comparison        
    CAP_holder_width = 10
    CAP_holder = DataFrame(bias=[], C1=[], C2=[], C3=[], C4=[], C5=[], C6=[], C7=[], C8=[], C9=[], C10=[])
  end 

  if SIM_list==Nothing
    if physical_model_name==Nothing
      physical_model_name = "cosi"
    end
    SIM_list = construct_SIM_list(TC=TC, pO2=pO2, bias=bias, data_set=data_set, simulations=["EIS"], 
                                    physical_model_name=physical_model_name)
  end
    
  for SIM in SIM_list, lambda_item in lambda
    # to add .... , tau_min_fac=tau_min_fac, tau_max_fac=tau_max_fac, tau_range_fac=tau_range_fac
    DRT_control = DRT_control_struct(lambda_item, tau_min_fac, tau_max_fac, tau_range_fac, peak_merge_tol, Nothing)
    SIM.DRT_control = DRT_control
    SIM.DRT_draw_semicircles=draw_semicircles
    SIM.DRT_backward_check=backward_check
    SIM.plot_option=plot_option
    SIM.f_range=f_range
    #@show SIM
    
    if mode=="EEC"
      EIS_df = EIS_get_RC_CPE_elements(R1, C1, R2, C2, alpha, R_ohm, f_range=f_range)
      plot_bool && typical_plot_sim(SIM, EIS_df, "! EEC ($R1, $C1) ($R2, $C2, $alpha)")
    elseif mode=="sim"
      EIS_df = ysz_fitting.simple_run([SIM], pyplot=(plot_bool ? 1 : 0), 
        prms_names=prms_names, 
        prms_values=prms_values, use_experiment=false)
    elseif mode=="exp"
      #@show SIM.checknodes
      EIS_df = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
      if SIM.data_set[end-1:end]==".z"
        plot_bool && typical_plot_exp(SIM, EIS_df, "!$(SIM.data_set)") 
      else
        plot_bool && typical_plot_exp(SIM, EIS_df) 
      end
    end

    if CAP_comparison
      DRT_actual = get_DRT(EIS_df, DRT_control)

      current_line = Array{Float32}(undef, CAP_holder_width+1)
      current_line[1] = SIM.bias
      for i in 1:CAP_holder_width
        if i <= length(DRT_actual.peaks_df.C)
          current_line[i+1] = DRT_actual.peaks_df.C[i]
        else
          current_line[i+1] = -1
        end        
      end
      push!(CAP_holder, current_line)
    end
    
    if backward_check
      DRT_actual = get_DRT(EIS_df, DRT_control)
      plot_bool && println("Fitness error = ",fitnessFunction(EIS_simulation(), DRT_actual.EIS_df, EIS_df))
      plot_bool && typical_plot_sim(EIS_simulation(800, 80, 0.0, use_DRT=false, DRT_control=DRT_control, plot_option=plot_option)..., DRT_actual.EIS_df, "! DRT_backward_check")
    end
  end
  
  if CAP_comparison    
    figure(CAP_plot_num)
    title("Capacitance VS Bias")
    xlabel("\$\\eta\$ [V]")
    ylabel("C [F]")
    for i in 1:CAP_holder_width
      current_C_range = []
      current_bias_range = []
      for j in 1:length(CAP_holder[!, Symbol("C$(i)")])
        CC = CAP_holder[!, Symbol("C$(i)")][j]
        if CC > -0.5
          append!(current_C_range, CC)
          append!(current_bias_range, CAP_holder[!, :bias][j])
        end
      end
      
      plot(current_bias_range, current_C_range, "-x")
    end 
    
    if CAP_plot_CAP_CV
      TC    = []
      pO2   = []   
      bias  = []
      for SIM in SIM_list
        push!(TC, SIM.TC)
        push!(pO2, SIM.pO2)
        push!(bias, SIM.bias)
      end
      
      prms_names_CAP = (prms_names..., CAP_bottleneck_prm)
      prms_values_CAP = (prms_values..., 1)
      ysz_fitting.simple_run(CAP_simulation(TC, pO2, analytical=false, voltrate=1.0e-5, upp_bound=maximum(bias), low_bound=minimum(bias)), 
                            pyplot=1, prms_names=prms_names_CAP, f_range=f_range, prms_values=prms_values_CAP, use_experiment=false)
      CAP_plot_num = gcf().number
    end

#       if CAP_plot_analytical
#         TODO !!!
#       end
  end
  
  return
end





###########################################################
###########################################################
#### Working space ########################################
###########################################################
###########################################################

function plot_R_ohm_dependence(;
        DD_list=[1].*1.0e-11, 
        nu_list=collect(0.1 : 0.1 : 0.9),
        TC=800)
  
  R_ohm_all = []  
  for DD in DD_list
    for nu in nu_list
      sigma, R_ohm = ysz_experiments.run_new(;
                pO2=1.0, T=TCtoT(TC), data_set="OLD_MONO_100",
                prms_names_in=["DD","nu", "weird_DD"],
                prms_values_in=(DD, nu, false),
                  
                physical_model_name="ysz_model_GAS_LoMA_Temperature",
                conductivity_fitting=true
                ) 
       push!(R_ohm_all, R_ohm)
    end
  end
  
  if length(nu_list) == 1
    plot(DD_list, R_ohm_all)
  else
    plot(nu_list, R_ohm_all)
  end
  return R_ohm_all
end


function model_R1_from_TC_sim(x, TC; print_DD=false)
  DD =  1.0e-13*x[1]*exp((e0*x[2])/(kB*TCtoT(TC)))

  if print_DD
    @show DD
  end
  sigma, R_ohm = ysz_experiments.run_new(;physical_model_name="",
              pO2=1.0, T=TCtoT(TC),
              prms_names_in=["DD","weird_DD_bool"],
              prms_values_in=(DD, false),
              
              conductivity_fitting=true
              )
  #@show x, R_ohm
  return R_ohm
end

function find_TC_fit_sim(TC_data, R1_data; initial_guess=[1.0e+7, -0.5])

  function to_optimize(x)
    error = 0
    for (i, TC) in enumerate(TC_data)
      error += (R1_data[i] - model_R1_from_TC_sim(x, TC))^2
    end
    #@show x, error
    #@show sqrt(error)/length(TC_data)
    return sqrt(error)/length(TC_data)
  end
  
  fit_O = optimize(to_optimize, initial_guess)
  
  R1_fitted_values = []
  TC_plot_range = 700 : 10 : 850
  for TC in TC_plot_range
    append!(R1_fitted_values, model_R1_from_TC_sim(fit_O.minimizer, TC))
  end
  
  title("Fitting of Diffusion Coefficient ... \$ D = A \\ \\mathrm{exp}(-\\frac{E}{k_\\mathrm{B} T})  \$")
  xlabel("TC")
  ylabel("\$R_{\\Omega}\$")
  plot(TC_data, R1_data, label="data: MONO_110, bias=0", "x")
  plot(TC_plot_range, R1_fitted_values, label="fit")
  legend(loc="best")
  return fit_O
end


###########################################################
###########################################################
#### Levenberg-Marquardt ##################################
###########################################################
###########################################################

function get_fitting_initial_condition()
  # fitted things for 850, 100, 0, MONO
  #x0 = (0.0, 0.0, 0.0, 20.01931310942335, 19.13133370088284, 20.30899006042895, 1, 21, -0.10582329539200844, -0.5612452246493124, -0.004425970436021694, 9.0e-13, 0.85, true, 0.25, 0.5)
  x0 = (0.0, 0.0, 0.0, 20.009046378779413, 19.133842639499733, 20.243362989853896, 1, 21, -0.13315570684220343, -0.5493240083421153, -0.01362036864906086, 9.0e-13, 0.85, true, 0.25, 0.5)
  #x0 =  (0.0, 0.0, 0.0, 20.304854396292857, 19.100695998282248, 20.203452081086848, 29.786448291451634, 22.50852142987179, -0.09991041679025696, -0.47107633036111757, -0.01598123040079823, 9.0e-13, 0.85, true, 0.25, 0.5)
  #x0 =  (0.0, 0.0, 0.0, 20.304854396292857, 19.100695998282248, 20.303452081086848, 25.786448291451634, 22.50852142987179, -0.09991041679025696, -0.47107633036111757, -0.01598123040079823, 9.0e-13, 0.85, true, 0.25, 0.5)
  
  x0 = (0.0, 0.0, 0.0, 19.85559447494341, 20.057302906172676, 19.30709527857658, 5.129458479008138, 25.342869613876285, -0.45254880302146255, -0.6377067007673729, -0.940802850242693, 9.0e-13, 0.85, true, 0.25, 0.5) # Optim ||  EIS: err =0.00442703203277445 (850, 100, 0.0)
  
  #x0 = (0.0, 0.0, 0.0, 19.907580949649912, 20.028132917464422, 19.403302626940977, 6.0097725379445714, 25.284531020122714, -0.46362112792111193, -0.6359175506862039, -0.9403837928159715, 9.0e-13, 0.85, true, 0.25, 0.5) # Optim ||  EIS: err =0.00442703203277445 (850, 100, 0.0)
  
  #x0 = (0.0, 0.0, 0.0, 19.7936, 20.2862, 18.7722, 5.0, 25.0, -0.508173, -0.648986, -0.795695, 9.0e-13, 0.85, true, 0.25, 0.5)
  
  x0 = (1.0, 0.0, 1.0, 20.85559447494341, 20.057302906172676, 20.30709527857658, 20.129458479008138, 20.342869613876285, -0.25254880302146255, -0.1377067007673729, -0.140802850242693, 9.0e-13, 0.85, true, 0.25, 0.5) # Optim ||  EIS: err =0.00442703203277445 (850, 100, 0.0)
  
  x0 = (1.0, 0.0, 1.0, 24.68638793247588, 16.929296123207436, 22.170921800982974, 18.553478252136934, 19.897488751668668, -0.5769408452228384, 0.12810885586649326, -0.33908733417798587, 9.0e-13, 0.85, true, 0.25, 0.5) # err = 0.03701 .. cosi pro (850, [60, 80, 100], 0.0) MONO
  
  
  x0 = (1.0, 0.0, 1.0, 23.9595, 23.4001, 21.1791, 5.0, 25.0, -0.78899, -0.522923, 0.591665, 9.0e-13, 0.85, true, 0.25, 0.5) # err =  0.025134653 .. cosi pro (850, [60, 80, 100], 0.0) MONO
  
  x0 = (1.0, 0.0, 1.0, 26.429724575986242, 21.599841316855407, 21.084721296249796, -6.663086428479204, 23.342499176410396, -2.123456036083893, -0.6371129649203712, 0.6658969910747834, 9.0e-13, 0.85, true, 0.25, 0.5) # err =  0.013466 .. cosi pro (850, [60, 80, 100], 0.0) MONO
  
  x0 = (1.0, 0.0, 1.0, 20.93358829626075, 21.402618660603686, 21.136328351766778, 5.47009852248677, 23.37319636187158, -0.4140996223040165, -0.6454107557668001, 0.6964208008004982, 9.0e-13, 0.85, true, 0.25, 0.5) # intermediate ... to delete !
  
  x0 = (1.0, 0.0, 1.0, 20.93358829626075, 20.402618660603686, 20.136328351766778, 0.0, 0.0, -0.0140996223040165, -0.06454107557668001, 0.06964208008004982, 12.3e-11, 0.85, true, 0.25, 0.5) # intermediate ... to delete !
  
  
  x0 = (1.0, 0.0, 1.0, 24.352, 25.9, 22.1599, 0.0, 0.0, 0.271907, -0.197024, -0.214102, 1.23e-10, 0.85, true, 0.25, 0.5)
  return x0
end

mutable struct SIM_fitting_struct
  TC
  pO2
  bias
  data_set
  simulations
  fitness_factors
  physical_model_name
  #
  SIM_list
  #
  prms_names::Array{String}
  x0
  mask
  lower_bounds
  upper_bounds
  prms_values
  #
  bboptimize_bool::Bool  #which optimizer is called ..... TODO!!! maybe change to string?
  iteration_count::Int64
  #
  print_to_file::Bool
  save_dir
  file_name
  
  SIM_fitting_struct() = new()
end

function build_SIM_fitting(;TC=850, pO2=80, bias=0.0, data_set="MONO_110", simulations=["EIS"], fitness_factors=[1.0],
                      physical_model_name="ysz_model_GAS_LoMA_Temperature",
                      #                      
                      prms_names=["kappaA", "kappaR", "kappaO", 
                                  "rA", "rR", "rO",         "rB", "rC",     
                                  "DGA", "DGR", "DGO",     
                                  "DD", "nu", "separate_vacancy",       "sites_Om0", "sites_Om1"  ],
                      x0 = get_fitting_initial_condition(),            
                      mask          =(0, 0, 0,
                                      1, 1, 1,        0, 0,
                                      1, 1, 1,
                                      0, 0,     0,        0, 0    ),
                      lower_bounds=(0.0, 0.0, 0.0,
                                    15.5, 15.9, 15.7,        5, 5,       
                                    -0.8, -0.8, -0.8,     
                                    [1]*1.0e-13, 0.01, true,       -Inf, -Inf    ),
                      upper_bounds=(0.0, 0.0, 0.0,
                                    25.5, 25.9, 25.7,        25, 25,        
                                    0.8, 0.8, 0.8,     
                                    [1]*1.0e-8, 0.99, true,       Inf, Inf    ),
                      #
                      bboptimize_bool=false, iteration_count=100,   
                      #
                      print_to_file=true, save_dir="../data/EEC/temp/log/", file_name="default.txt",
                      )
  SIM_fitting = SIM_fitting_struct()
  
  SIM_fitting.TC = TC
  SIM_fitting.pO2 = pO2
  SIM_fitting.bias = bias
  SIM_fitting.data_set = data_set
  SIM_fitting.simulations = simulations
  SIM_fitting.fitness_factors = fitness_factors
  SIM_fitting.physical_model_name = physical_model_name
  #
  SIM_fitting.prms_names = prms_names
  #@show typeof(x0)
  SIM_fitting.x0 = x0
  SIM_fitting.mask = mask
  SIM_fitting.lower_bounds = lower_bounds
  SIM_fitting.upper_bounds = upper_bounds
  #
  SIM_fitting.bboptimize_bool = false
  SIM_fitting.iteration_count = iteration_count
  #
  SIM_fitting.print_to_file = print_to_file
  SIM_fitting.save_dir = save_dir
  SIM_fitting.file_name = file_name
  
  SIM_fitting.SIM_list = get_SIM_list_rectangle(TC, pO2, bias, data_set, simulations, fitness_factors, physical_model_name)

  return SIM_fitting
end

function run_SIM_fitting_script_wrap(
                    TC_string, 
                    pO2_string,
                    bias_string,
                    data_set,
                    simulations_string,                    
                    fitness_factors_string,
                    physical_model_name,
                    #                    
                    prms_names_string,
                    x0_string,            
                    mask_string,
                    lower_bounds_string,
                    upper_bounds_string,
                    #
                    bboptimize_bool_string, 
                    iteration_count_string,   
                    #
                    print_to_file_string, 
                    save_dir, 
                    file_name,
                    #
                    #
                    pyplot_string,  
                    plot_each_x_th_string,
                    print_only_result_string,
                    )
  
  SIM_fitting = SIM_fitting_struct()
  #
  TC = eval(Meta.parse(TC_string))
  pO2 = eval(Meta.parse(pO2_string))
  bias = eval(Meta.parse(bias_string))
  fitness_factors = eval(Meta.parse(fitness_factors_string))
  simulations = eval(Meta.parse(simulations_string))
  #
  SIM_fitting.TC = TC
  SIM_fitting.pO2 = pO2
  SIM_fitting.bias = bias
  SIM_fitting.data_set = data_set
  SIM_fitting.simulations = simulations
  SIM_fitting.fitness_factors = fitness_factors
  SIM_fitting.physical_model_name = physical_model_name
  #
  SIM_fitting.SIM_list = get_SIM_list_rectangle(TC, pO2, bias, data_set, simulations, fitness_factors, physical_model_name)
  #
  SIM_fitting.prms_names = eval(Meta.parse(prms_names_string))
  SIM_fitting.x0 = eval(Meta.parse(x0_string))
  SIM_fitting.mask = eval(Meta.parse(mask_string))
  SIM_fitting.lower_bounds = eval(Meta.parse(lower_bounds_string))
  SIM_fitting.upper_bounds = eval(Meta.parse(upper_bounds_string))
  #
  SIM_fitting.bboptimize_bool = eval(Meta.parse(bboptimize_bool_string))
  SIM_fitting.iteration_count = eval(Meta.parse(iteration_count_string))
  #
  SIM_fitting.print_to_file = eval(Meta.parse(print_to_file_string))
  SIM_fitting.save_dir = save_dir
  SIM_fitting.file_name = file_name
  
  
  
  printfields(SIM_fitting)
  return run_SIM_fitting(SIM_fitting, 
                  pyplot=eval(Meta.parse(pyplot_string)),
                  plot_each_x_th=eval(Meta.parse(plot_each_x_th_string)),
                  print_only_result=eval(Meta.parse(print_only_result_string))
                  )
end

# function SIM_fitting_save_to_file(SIM_fitting::SIM_fitting_struct, )
#   
# end


function SIM_fitting_show_EIS(SIM_fitting::SIM_fitting_struct, pyplot=1, use_experiment=true, use_fitted_values=false)
  simple_run(
    SIM_fitting.SIM_list, 
    prms_names=SIM_fitting.prms_names,
    prms_values= (use_fitted_values ? SIM_fitting.prms_values : SIM_fitting.x0),
    use_experiment=use_experiment,
    pyplot=pyplot
    )
end

function run_SIM_fitting_example()
  SIM_fitting = build_SIM_fitting()
  
  fit = run_SIM_fitting(SIM_fitting, pyplot=false, print_only_result=true)
  return fit
end


function run_SIM_fitting(SIM_fitting::SIM_fitting_struct;
                      pyplot=false,  plot_each_x_th=50,
                      print_only_result=false,
                      )
  function plot_error_projection(prms_values, prms_names, error)
    fig_num = 333
    prms_length = length(prms_names)
    plot_edge = ceil(sqrt(prms_length))
    
    if projection_plot_maximum < 0 && error < 20
      projection_plot_maximum = error
    end
    
    if error < projection_plot_minimum
      projection_plot_minimum = error
    end
    
    figure(fig_num)
    for i in 1:prms_length
      subplot(plot_edge, plot_edge, i)
      xlabel(prms_names[i])
      ylabel("err")
      ylim(projection_plot_minimum, projection_plot_maximum)
      plot(prms_values[i], error, "o")
    end
    
    return
  end
  
#   function prepare_prms(mask, x0, x)
#       prms = []
#       xi = 1
#       for i in collect(1 : 1 :length(mask))
#           if convert(Bool,mask[i])
#               append!(prms, x[xi])
#               xi += 1
#           else
#               append!(prms, x0[i])
#           end
#       end
#       return Tuple(prms)
#   end
  
  function to_optimize(x)
      #@show x
      if !(check_x_in(x, lowM, uppM))
        print_only_result || (output_buffer *= "    OUT OF THE BOUNDS   \n")
        print_only_result || (print_and_delete_buffer(output_buffer))
        return 10000
      end
      iteration_counter += 1
      
      prms_values = prepare_prms(SIM_fitting.mask, SIM_fitting.x0, x)
#       print(" >> SIM_fitting.mask = ",SIM_fitting.mask)
#       print(" || prms = ",prms_values)

      print_only_result || (output_buffer *= "> prms=$(prms_values)\n")
      err = 0.0
      SIM_err = 0.0
      err_string = ""

      
      for (i, SIM) in enumerate(SIM_fitting.SIM_list)        
        try
          if check_equal_size(SIM_fitting.prms_names, prms_values)  
            SIM_sim = apply_checknodes(SIM, 
                      typical_run_simulation(SIM, SIM_fitting.prms_names, prms_values, pyplot ? 1 : 0),
                      SIM.checknodes)
          end    
                
          if pyplot && mod(iteration_counter, plot_each_x_th) == 1
            typical_plot_sim(SIM, SIM_sim)
          end
          
          SIM_err = SIM.fitness_factor * fitnessFunction(SIM, SIM_sim, SIM_exp[i])
          err += SIM_err
        catch e
          if e isa InterruptException
            rethrow(e)
          else
            print(e)
            if norm(x0M.-x) < 0.01
              err += 10000
            else
              err += 10000*norm(x0M.-x) 
            end
          end
        end 
        
        err_string = " "*string(SIM)
        print_only_result || (output_buffer *= "$(err_string)=$(SIM_err) || ")
      end
      print_only_result || (output_buffer *= "err=$(err)\n\n")
      
      
      if pyplot
        plot_error_projection(
          take_only_masked(SIM_fitting.mask, prms_values),
          take_only_masked(SIM_fitting.mask, SIM_fitting.prms_names), 
          err)
      end
      
      print_only_result || (print_and_delete_buffer(output_buffer))
      
#       println(" || ",err_string,": err =", err)
      println(iteration_counter, " >> x =", round.(x, digits=3)," : err =", err)
      return err
  end
  
  function print_and_delete_buffer(buffer="string")
    if SIM_fitting.print_to_file
      open(SIM_fitting.save_dir*SIM_fitting.file_name, "a") do f
        write(f, buffer)
      end
    else
      print(buffer)
    end
  end

  
  println(" -------- run_SIM_fitting --- started! -------- ")
  #return 0
    
  iteration_counter = 0
  
  if pyplot
      PyPlot.close()
  end
  
  SIM_exp = Array{DataFrame}(undef, length(SIM_fitting.SIM_list))
  for (i, SIM) in enumerate(SIM_fitting.SIM_list)
    SIM.plot_legend=false
    SIM_exp[i] = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
    if (pyplot > 0)
      typical_plot_exp(SIM, SIM_exp[i])
    end
  end

  output_buffer = ""
  if SIM_fitting.print_to_file
    mkpath(SIM_fitting.save_dir)
  end
  #######################################
  #######################################
  #######################################
  #######################################
  # control panel
  #SIM_fitting.x0 = zeros(2)
  #optimize(rosenbrock, SIM_fitting.x0, LevenbergMarquardt())
  

  ######################################            

  if length(SIM_fitting.x0) != length(SIM_fitting.prms_names)
    output_buffer *= ("ERROR: length(SIM_fitting.x0) != length(SIM_fitting.prms_names)")
    print_and_delete_buffer(output_buffer)
    return throw(Exception)
  end

  
  #######################################
  #######################################
  #######################################

  #######################################
  #######################################
  #######################################
  #######################################
  
  projection_plot_maximum = -1
  projection_plot_minimum = Inf
  
  x0M = zeros(0)
  lowM = zeros(0)
  uppM = zeros(0)
  for i in collect(1 : 1 : length(SIM_fitting.mask))
      if convert(Bool,SIM_fitting.mask[i])
          append!(x0M, SIM_fitting.x0[i])
          append!(lowM, SIM_fitting.lower_bounds[i])
          append!(uppM, SIM_fitting.upper_bounds[i])
      end
  end
  
  function get_SearchRange(lowM, uppM)
    output = Array{Tuple{Float64, Float64}}(undef, 0)
    for i in 1:length(lowM)
      append!(output, [(lowM[i], uppM[i])])
    end
    return output
  end
  
  # for BlackBoxOptim
  SearchRange = get_SearchRange(lowM, uppM)
  #@show SearchRange
  
  if SIM_fitting.bboptimize_bool
    fit= bboptimize(to_optimize, SearchRange = SearchRange)
  else
#     fit = optimize(
#       to_optimize,
#       x0M,
#       #rosenbrock2d,
#       #[0., 0.],
#       #NelderMead(),
#       ParticleSwarm(lower=lowM, upper=uppM, n_particles = 10) 
#       )
      
    fit = optimize(to_optimize, x0M, NelderMead(),
                     Optim.Options(
                      iterations = SIM_fitting.iteration_count,
                      
                      )
                     )
     
#    fit = optimize(to_optimize, lowM, uppM, x0M, Fminbox(NelderMead()), 
#                     Optim.Options(
#                      iterations = 10,
#                      
#                       )
#                     )
  end
  
  SIM_fitting.prms_values = prepare_prms(SIM_fitting.mask, SIM_fitting.x0, fit.minimizer)
  
  if print_only_result
    output_buffer *= "\n"
    output_buffer *= "SIM_list = $(string.(SIM_fitting.SIM_list))\n"
    output_buffer *= "data_set = $(string.(SIM_fitting.data_set))\n"
    output_buffer *= "prms_names = $(SIM_fitting.prms_names)\n"
    output_buffer *= "x0 = $(SIM_fitting.x0)\n"
    output_buffer *= "lower_bounds = $(SIM_fitting.lower_bounds)\n"
    output_buffer *= "upper_bounds = $(SIM_fitting.upper_bounds)\n"
    output_buffer *= "prms_values = $(SIM_fitting.prms_values)\n"
    output_buffer *= "mask = $(SIM_fitting.mask)\n"
    output_buffer *= "err=$(fit.minimum)\n"
    print_and_delete_buffer(output_buffer)
  end
  
  #println(optimize(to_optimize, x0M, lower=lowM, upper=uppM, Δ=1000000, f_tol=1.0e-14, g_tol=1.0e-14, LevenbergMarquardt()))
  
  ####optimize(to_optimize, x0M, Dogleg())
  return fit  
end


function partial_meta_run_SIM_fitting(;SIM_fitting::SIM_fitting_struct=build_SIM_fitting(),
                              pyplot = false,
                              plot_each_x_th = 50,
                              print_only_result = true,
                              direct_bool = true,
                              bash_command = "julia",
                              run_file_name = "../snehurka/run_ysz_fitting_SIM_fitting.jl"      
                              )
  if direct_bool
    return run_SIM_fitting(SIM_fitting,
                      pyplot=pyplot,
                      plot_each_x_th=plot_each_x_th,
                      print_only_result=print_only_result
                      )
  else
    
    return run(`
                $bash_command  
                $run_file_name 
                
                $(string(SIM_fitting.TC)) 
                $(string(SIM_fitting.pO2))
                $(string(SIM_fitting.bias))
                $(string(SIM_fitting.data_set))
                $(string(SIM_fitting.simulations))
                $(string(SIM_fitting.fitness_factors))
                $(string(SIM_fitting.physical_model_name))
                
                                    
                $(string(SIM_fitting.prms_names))
                $(string(SIM_fitting.x0))
                $(string(SIM_fitting.mask))
                $(string(SIM_fitting.lower_bounds))
                $(string(SIM_fitting.upper_bounds))
                
                $(string(SIM_fitting.bboptimize_bool))
                $(string(SIM_fitting.iteration_count))
                
                $(string(SIM_fitting.print_to_file))
                $(string(SIM_fitting.save_dir))
                $(string(SIM_fitting.file_name))
                
                
                $(string(pyplot))
                $(string(plot_each_x_th))
                $(string(print_only_result))
    `)
  end
end


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################



function slurm_evaluate_results(;print_bool=false, show_x0 = false, working_dir="../snehurka/", search_last_lines=400)

  function SER_get_value(line, name; set_standard_failed=true)
    error_split = split(line, '=')
    if occursin(name, error_split[1])
      num = eval(Meta.parse(
            error_split[2]
            )
          )
      if typeof(num)!=Nothing
      	return num
      else
        return ""
      end
    else
      standard_failed = true
      return ""
    end
  end


  dir_items = cd(readdir, working_dir)

  standard_failed = false
  non_standard_success = false
  error_and_prms_values_dataframe = DataFrame(error = [], prms_values = [], x0 = [], file_name = [])
  for item in dir_items
    if (length(item) >= 5) && (item[1:5]=="slurm")
      standard_failed = false
      non_standard_success = false

      buffer = readlines(working_dir*item)

      if length(buffer) > 3
        error = SER_get_value(buffer[end-1], "err")
        prms_values = SER_get_value(buffer[end-3], "prms_values")
        x0 = SER_get_value(buffer[end-6], "x0")


        if standard_failed
          mask = []
          x0 = []
          x_values = []

          for line in buffer
            pre_mask = SER_get_value(line, "mask", set_standard_failed=false)
            pre_mask=="" ? false : mask = pre_mask
            pre_x0 = SER_get_value(line, "x0", set_standard_failed=false)
            pre_x0=="" ? false : x0 = pre_x0

            if (length(mask) > 0) & (length(x0) > 0)
              break
            end
          end
          if (length(mask) < 1) || (length(x0) < 1)
            continue
          end

          error = Inf
          actual_error = 0
          starting_line_number = max(length(buffer) - search_last_lines, 1)
          for line in buffer[starting_line_number : end]
            aux = split(line, '=')

            first_raw = split(aux[1], ">>")
            if length(first_raw) > 1
              x_str = first_raw[2]
            else
              continue
            end
            if x_str == " x "
              try
                x_values = eval(Meta.parse(split(aux[2], ':')[1]))
                actual_error = eval(Meta.parse(aux[3]))
                if length(x_values) < 1 || length(actual_error) < 1
                  continue
                else
                  non_standard_success = true
                end
              catch
                continue
              end
            else
              continue
            end
            if actual_error < error
              error=actual_error
              prms_values = prepare_prms(mask, x0, x_values)
            end
          end

        end

        if !(standard_failed) || (non_standard_success)
          push!(
            error_and_prms_values_dataframe,
            (error, prms_values, x0, item)
          )
        end
      end
    end
  end

  sort!(error_and_prms_values_dataframe, :error)
  if print_bool
    if !show_x0
      @show error_and_prms_values_dataframe[:, [:error, :prms_values]]
    else
      @show error_and_prms_values_dataframe[:, [:error, :x0]]
    end
  end
  return error_and_prms_values_dataframe
end



function SIM_fitting_evaluate(SIM_fitting, fitted_values)
  simple_run(
    SIM_fitting.SIM_list,
    prms_names=SIM_fitting.prms_names,
    prms_values= fitted_values,
    use_experiment=true,
    pyplot=1
    )
end

function SIM_fitting_x0_test(SIM_fitting, x0; pyplot=false, print_only_result=false, plot_each_x_th=45)
  SIM_fitting.x0 = x0
  return run_SIM_fitting(SIM_fitting, pyplot=pyplot, print_only_result=print_only_result, plot_each_x_th=plot_each_x_th)
end


function aux_save_SIM_fitting(res::DataFrame, name)
  file_path = "aux_slurm_data/"*name*"/"
  run(`mkdir -p $(file_path)`)
  
  dir_items = cd(readdir, working_dir)  
  #
  # TODO!!
  #
end



function plot_temp_parameters(;prms_names, prms_values, TC_list = [700, 750, 800, 850],
                              show_reactions=["A","R","O"], 
                              show_parameters=["r", "DG", "nu", "CO", "COmm"],
                              label="_set_1",
                              save_file=false, save_dir="../data/temp_prms/", file_name="temp_prms_$(label)",                              
                              fig_num = 135,
                              f_b_rates=true,
                              line_color_idx = 1,                              
                              )
  line_styles = ["-", "--", ":"]
  line_colors = ["r", "g", "b"]
    
    
  function gv(prm_name; throw_error=false)
    idx_tuple = findall(x->x==prm_name, prms_names)
    if length(idx_tuple) == 0
      if throw_error
        println("ERROR: parameter $(prm_name) not found!")
        return throw(Exception)
      else
        return Nothing
      end
    else
      return prms_values[idx_tuple[1]]
    end
  end
  
  function get_X_list(prm_name, TC_list)
    X_B = gv(prm_name*"_B")
    X_C = gv(prm_name*"_C")
    
    #
    # TODO !!!!!!! >> use appropriate model and its interpretation of parameters
    #
    if X_B == Nothing || X_C == Nothing
      return [gv(prm_name*"_$(TC)", throw_error=true) for TC in TC_list]
    else
      if (length(prm_name) >= 3) && (prm_name[2:3]==".r")        
        return [10.0^((X_B*((TC-700)/50.) + X_C) + (gv(prm_name[1]*".S"))) for TC in TC_list]
      else
        return [X_B*((TC-700)/50.) + X_C for TC in TC_list]
      end
    end
  end
  
  function plot_X(prm_name, special_y_label=Nothing, special_legend=Nothing; line_style=Nothing)
    xlabel("TC (°C)")
    special_y_label == Nothing ? ylabel(prm_name) : ylabel(special_y_label)
    
    #subplots_adjust(hspace = 0.5)
    
    X_list = get_X_list(prm_name, TC_list)
    if save_file
      TC_prms_df[!, Symbol(prm_name)] = X_list
    end
    
    if (length(prm_name) >= 3) && (prm_name[2:3]==".r")
      yscale("log")
    end
    
    plot(TC_list, X_list, 
      label= (special_legend==Nothing ? prm_name*"_"*label : 
        (special_legend=="!none" ? "" : special_legend)
      ),
      (line_style == Nothing ? "-" : line_style)
    )  
    legend(loc="best")
    grid(true)
    return X_list
  end
  
  
  
  
  include("../src/models/ysz_model_GAS_LoMA_Temperature.jl")
  physical_model_name="ysz_model_GAS_LoMA_Temperature"
  
  #model_symbol = eval(Symbol("ysz_experiments."*physical_model_name))
  AreaEllyt = 0.00011309724 * 0.7
  parameters=ysz_experiments.ysz_model_GAS_LoMA_Temperature.YSZParameters()
  
  # formal experimental setting for pO2
  parameters.pO2 = 0.20
  

  


  
  
  figure(fig_num)
  suptitle("Temperature dependent parameters")
  
  if save_file 
   TC_prms_df = DataFrame(TC = TC_list)
  end
    
  r_idx = 1
  DG_idx = 2
  
  T_list = TCtoT.(TC_list)

  reactions_lists=[[],[]]
  
  for (i, reaction_name) in enumerate(show_reactions)
    for (j, attribute_name) in enumerate(["r", "DG"])
      prm_name = reaction_name*"."*attribute_name      
      
      subplot(2, 2, j)
      push!(reactions_lists[j],
        plot_X(
            prm_name, 
            (attribute_name=="r" ? "r [m/s^2]" : "DG [eV]"), 
            (line_color_idx == 1 ? reaction_name : "!none"), 
            line_style=line_colors[line_color_idx]*line_styles[i]
        )    
      )
    end  
  end
  
  subplot(3, 3, 6 + 1)
  nu_list = plot_X("nu", line_style=line_colors[line_color_idx])
  
  subplot(3, 3, 6 + 2)
  CO_list = plot_X("CO", line_style=line_colors[line_color_idx])
  
  subplot(3, 3, 6 + 3)
  COmm_list = plot_X("COmm", line_style=line_colors[line_color_idx])
  
  # plot forward and backward rates
  if f_b_rates
    figure(fig_num + 1)
    suptitle("Derived quantities")
    
    

    subplot(2,2,2)
    plot([], line_colors[line_color_idx], label=label)
    legend(loc="best")
    
    
    
    r_f_lists = []
    r_b_lists = []
    for (i, reaction_name) in enumerate(show_reactions)
      # forward
      push!(r_f_lists, 
        (reactions_lists[r_idx][i]/(10.0^gv(reaction_name*".S"))) .* exp.(- (10.0^gv(reaction_name*".S"))* reactions_lists[DG_idx][i] .* eV ./ (2*kB*T_list))
      )
      subplot(2,2,1)
      
# #       cm=get_cmap(:tab20)
# #       cycler = pyimport("cycler")
# #       PyPlot.rc("axes",prop_cycle=cycler.cycler(color=[cm(t/19) for t in 0:19]))

      xlabel("TC (°C)")
      ylabel("r_forward")
      yscale("log")
      plot( 
        TC_list, r_f_lists[i], 
        line_colors[line_color_idx]*line_styles[i], 
        #label=reaction_name
        label=(line_color_idx == 1 ? reaction_name : "")
        )
      legend(loc="best")
      if save_file
        TC_prms_df[!, Symbol("r_f_"*reaction_name)] = r_f_lists[end]
      end
      
      # backward
      push!(r_b_lists, 
        (reactions_lists[r_idx][i]/(10.0^gv(reaction_name*".S"))) .* exp.( (10.0^gv(reaction_name*".S"))*reactions_lists[DG_idx][i] .* eV ./ (2*kB*T_list))
      )
      subplot(2,2,2)
      xlabel("TC (°C)")
      ylabel("r_backward")
      yscale("log")
      plot( TC_list, r_b_lists[i], line_colors[line_color_idx]*line_styles[i]
        #, label=reaction_name
        )
      #legend(false)
      #legend(loc="best")
      if save_file
        TC_prms_df[!, Symbol("r_b_"*reaction_name)] = r_b_lists[end]
      end
      
    end
  end
  
  function local_update_physical_parameters(TC)
    # update parameters for new TC
    parameters.T = TCtoT(TC)
    ysz_experiments.ysz_model_GAS_LoMA_Temperature.set_parameters!(parameters, prms_values, prms_names)
    return ysz_experiments.ysz_model_GAS_LoMA_Temperature.YSZParameters_update!(parameters)  
  end
  
  E_eq_list = []
  YSZ_Om0_eq_list = []
  surf_Om_eq_list = []
  surf_O_eq_list = []
  
  for TC in TC_list
    parameters = local_update_physical_parameters(TC)
    push!(E_eq_list, parameters.phi_eq*(1 + parameters.e_fac))
    push!(YSZ_Om0_eq_list, parameters.y0_eq*(1-parameters.nu)*parameters.m_par)
    push!(surf_Om_eq_list, parameters.yAs_eq*parameters.COmm)
    push!(surf_O_eq_list, parameters.yOs_eq*parameters.CO)
  end
  
  
  function mymy_plot(X_list, yylabel, position)
    subplot(2,4, 4 + position)
    xlabel("TC (°C)")
    title(yylabel)
    plot(
      TC_list, X_list, 
      line_colors[line_color_idx],
      label=label
      )
    legend(loc="best")
    if save_file
      TC_prms_df[!, Symbol(yylabel)] = X_list
    end
  end
  
  mymy_plot(E_eq_list, "E_eq", 1)
  mymy_plot(YSZ_Om0_eq_list, "YSZ_Om0_eq", 2)
  mymy_plot(surf_Om_eq_list, "surf_Om_eq", 3)
  mymy_plot(surf_O_eq_list, "surf_O_eq", 4)
  
  
  if save_file
    mkpath(save_dir)
    CSV.write(save_dir*file_name*".csv", TC_prms_df)
  end
  
  return 
end



###########################################################
###########################################################
#### par_search_SIM_fitting ###############################
###########################################################
###########################################################

# function par_search




###########################################################
###########################################################
#### CLUSTER COMPUTATION ##################################
###########################################################
###########################################################




function meta_run_par_study(;only_return_SIM_fitting=false,
                            prms_lists,
                            name="default_name",
                            pyplot=false,
                            plot_each_x_th=20,
                            print_only_result=true,
                            SIM_fitting=Nothing,
                            scripted_tuple,
                            
                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = true,
  
                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!
                  
                            bash_command = "sbatch",
                            #bash_command = "echo",
                            #bash_command = "julia",
                            
                            #mode = "test_one_prms",
                            #mode = "only_print",
                            mode = "go",
                            
                            express3_bool = true
                            )  
                            
  function recursive_bash_call(output_prms_lists, active_idx)
    if active_idx > size(scripted_tuple,1)
      scripted_tuple_string = string(scripted_tuple)
      output_prms_lists_string = string(output_prms_lists)
      prms_names_string = string(SIM_fitting.prms_names)

      TC_string = string(SIM_fitting.TC)
      pO2_string = string(SIM_fitting.pO2)
      bias_string = string(SIM_fitting.bias) 
      data_set = SIM_fitting.data_set
      simulations_string = string(SIM_fitting.simulations)
      fitness_factors_string = string(SIM_fitting.fitness_factors)
      #physical_model_name = physical_model_name
      
      ### nasty merge ...... <<<<<<<<<<<<< TODO TODO TODO !!!!!!!!!!!!!!!! >>> make separate function for SIM_fitting
      if SIM_fitting_mode
        
        #@show output_prms_lists
        SIM_fitting.x0 = output_prms_lists
        
        partial_meta_run_SIM_fitting(SIM_fitting=SIM_fitting,
                              pyplot=pyplot,
                              plot_each_x_th=plot_each_x_th,
                              print_only_result =print_only_result,
                              direct_bool = direct_bool,
                              bash_command=bash_command,
                              run_file_name=run_file_name
                              )
        return
        
                                    
      else
        if direct_bool
          return run_par_study_script_wrap(
                      output_prms_lists_string ,
                      save_dir ,
                      name ,
                      scripted_tuple_string ,
                      prms_names_string ,
                      TC_string ,
                      pO2_string ,
                      bias_string ,
                      simulations_string ,
                      mode ,
                      physical_model_name        
          )
        else
          return run(`
                      $bash_command  
                      $run_file_name 
                      $output_prms_lists_string 
                      $save_dir 
                      $name 
                      $scripted_tuple_string 
                      $prms_names_string 
                      $TC_string 
                      $(pO2_string) 
                      $bias_string 
                      $simulations_string 
                      $mode 
                      $physical_model_name
          `)
        end
      end
    else
      if scripted_tuple[active_idx] == 1
        for i in prms_lists[active_idx]
          recursive_bash_call(push!(deepcopy(output_prms_lists),i), active_idx + 1)
        end
      else
        recursive_bash_call(push!(deepcopy(output_prms_lists), prms_lists[active_idx]), active_idx + 1)
      end
    end
  end
  
  function consistency_check()
    if  (size(SIM_fitting.prms_names,1) != size(prms_lists,1)) ||
        (size(SIM_fitting.prms_names,1) != size(SIM_fitting.mask,1)) || 
        (size(SIM_fitting.prms_names,1) != size(SIM_fitting.lower_bounds,1)) || 
        (size(SIM_fitting.prms_names,1) != size(SIM_fitting.upper_bounds,1)) || 
        (size(SIM_fitting.prms_names,1) != size(scripted_tuple,1))
        
      println("ERROR: shape mismatch: lengths of prms_lists, mask, lower_bounds, upper_bounds and scripted_tuple are NOT the same")
      return throw(Exception)
    end
    
    my_eps = 1.0e-5
    for (i, prms_lists_item) in enumerate(prms_lists)
      if length(prms_lists_item) < 1
        println("ERROR: empty list for \"$(SIM_fitting.prms_names[i])\" in prms_lists")
        return throw(Exception)
      end
      for (j, prm) in enumerate(prms_lists_item)
        if  (prm < (SIM_fitting.lower_bounds[i] - my_eps)) || 
            (prm > (SIM_fitting.upper_bounds[i] + my_eps))
          println("ERROR: value $(prm) of parameter \"$(SIM_fitting.prms_names[i])\" is NOT in bounds [$(SIM_fitting.lower_bounds[i]), $(SIM_fitting.upper_bounds[i])]!")
          return throw(Exception)        
        end
      end
    end        
    
    return true
  end


  
  ##### TODO!!!! tohle bych mohl udelat taky obecne, at skrz jeden skriptovaci soubor muze projit vse
  
  if SIM_fitting_mode
    if express3_bool
      run_file_name = "../snehurka/run_EX3_ysz_fitting_SIM_fitting.jl"
    else
      run_file_name = "../snehurka/run_ysz_fitting_SIM_fitting.jl"
    end
  else
    if express3_bool
      run_file_name = "../snehurka/run_EX3_ysz_fitting_par_study-prms-.jl"
    else
      run_file_name = "../snehurka/run_ysz_fitting_par_study-prms-.jl"
    end  
  end
  
  #######################################################
  ####################################################### 
  #######################################################
  #######################################################
  consistency_check()
  if  (mode != "test_one_prms" &&
      mode != "only_print" &&
      mode != "go")
    println("ERROR: unkonwn mode \"$(mode)\"")
    return throw(Exception)
  end
  
  if mode == "test_one_prms"
    prms_lists = [list[Int64(ceil(end/2.0))] for list in prms_lists]
  end
  
  # counter of output files
  nodes_count = 1
  per_node_count = 1
  for (i, byte) in enumerate(scripted_tuple)
    if byte == 1
      nodes_count *= size((prms_lists)[i],1)
    else
      per_node_count *= size((prms_lists)[i],1)
    end 
  end
  
  
  # saving metafile

  metafile_string = "-----------  METAFILE for par_study  -------------\n"
  metafile_string = string(metafile_string,"name=", name,"\n")
  metafile_string = string(metafile_string,"physical_model_name=", SIM_fitting.physical_model_name,"\n")
  prms_strings = [string(item) for item in prms_lists]
  for (i,str) in enumerate(prms_strings)
    metafile_string = string(metafile_string, SIM_fitting.prms_names[i],"_list=",str,"\n")
  end
  
  metafile_string = string(metafile_string,"#### nodes / pernode_count  =  ", nodes_count," / ",per_node_count," ( = ",per_node_count*3.0/60.0,"m)  #### \n")
  metafile_string = string(metafile_string,"prms_names=", SIM_fitting.prms_names,"\n")
  metafile_string = string(metafile_string,"scripted_tuple=", scripted_tuple,"\n")
  metafile_string = string(metafile_string,"TC=", SIM_fitting.TC,"\n")
  metafile_string = string(metafile_string,"pO2=", SIM_fitting.pO2,"\n")
  metafile_string = string(metafile_string,"bias=", SIM_fitting.bias,"\n")
  metafile_string = string(metafile_string,"data_set=", SIM_fitting.data_set,"\n")
  metafile_string = string(metafile_string,"simulations=", SIM_fitting.simulations,"\n")
  #metafile_string = string(metafile_string,"string(SIM_list)=", string(SIM_list),"\n") 
  metafile_string = string(metafile_string,"save_dir=", SIM_fitting.save_dir,"\n")
  
  metafile_string = string(metafile_string,"mode=", mode,"\n")
  metafile_string = string(metafile_string,"bash_command=", bash_command,"\n")
  metafile_string = string(metafile_string,"run_file_name=", run_file_name,"\n")
  metafile_string = string(metafile_string,"prms_lists=", prms_lists,"\n")
  #metafile_string = string(metafile_string,"SIM_list=", SIM_list,"\n")
  
  
  
  if mode == "only_print"
    println(metafile_string)
    return
  end
  
  meta_file_path = string(SIM_fitting.save_dir, name, "/")
  run(`mkdir -p $(meta_file_path)`)
  write("$(meta_file_path)__metafile_par_study.txt", metafile_string)
  
  println()
  println(metafile_string)
  
  #setting jobs
  recursive_bash_call([], 1)
    
  println("ok :) ")
  
  if SIM_fitting_mode
    return SIM_fitting
  end
end
        


end
