module ysz_experiments
### WHAT was done since last commit ###############
# * phi_out = (phi_prev + phi)/2.0
# * AreaEllyt added to EIS current computation
# 
# TODO !!!
# [x] ramp for EIS steady state
# [x] better initial conditions
###################################################






using Printf
using VoronoiFVM
using PyPlot
using DataFrames
using CSV
using LeastSquaresOptim
using LinearAlgebra
using SparseArrays
using Plots

##########################################
# internal import of YSZ repo ############
#model_label = "ysz_model_GAS_exp_ads"
#model_label = "ysz_model_GAS_LoMA"
#model_label = "ysz_model_GAS_LoMA_shared"
#model_label = "ysz_shn"
#model_label = "ysz_shn_overvoltage"
#model_label = "ysz_shn_overvoltage__TEMP"

#include("../src/models/$(model_label).jl")

include("../src/models/ysz_model_GAS_LoMA.jl")
include("../src/models/ysz_model_GAS_LoMA_shared.jl")
include("../src/models/ysz_model_GAS_LoMA_Temperature.jl")
include("../src/models/ysz_model_GAS_LoMA_generic.jl")
include("../src/models/YSZ_LSM_model_ARO.jl")


include("../prototypes/timedomain_impedance.jl")



# --------- end of YSZ import ---------- #
##########################################



function run_new(;physical_model_name,
                test=false, test_from_above=false, print_bool=false, debug_print_bool=false, out_df_bool=false,
                verbose=false, pyplot=false, pyplot_finall=false, save_files=false,
                width=0.0005, dx_exp=-9,
                pO2=1.0, T=1073, data_set=Nothing,
                prms_names_in=Nothing,
                prms_values_in=Nothing,
                #### EIS ####
                EIS_IS=false,  bias=0.0, f_range=(0.9, 1.0e+5, 1.1),
                #
                EIS_TDS=false, tref=0,
                #### CV #####
                voltammetry=false, voltrate=0.010, upp_bound_eta=1.0, low_bound_eta=-1.0, sample=30, checknodes=[], 
                #### CV(f) ##
                fast_CV_mode=false,
                #### CAP ####
                dlcap=false, dlcap_analytical=false, 
                #### STEP ####
                STEP=false, dt_fac=1.5, dt_start=1.0e-5, t_end=0.01,
                #### conductivity ####
                conductivity_fitting=false
                )

                
                
                
    # Geometry of the problem
    #AreaEllyt = 0.000201 * 0.6      # m^2   (geometrical area)*(1 - porosity)
    
    #porosity = 0.0
    porosity = 0.3
    if data_set == Nothing
      AreaEllyt = 1.0    # random value
    else
      if data_set == "OLD_MONO_100"
        AreaEllyt = 0.00011309724 * (1-porosity)       # m^2 (geometrical area)*(1 - porosity)
      elseif data_set[1:4] == "MONO"
        AreaEllyt = 0.00011309724 * (1-porosity)       # m^2 (geometrical area)*(1 - porosity)
      elseif data_set[1:4] == "POLY"
        AreaEllyt = 0.000201 * (1-porosity)            # m^2 (geometrical area)*(1 - porosity)
      else
        println("\nERROR: ysz_experiments.jl: data_set NOT RECOGNIZED !!!! \n")
        return throw(Exception)
      end
    end             
                
                
                
                
    ### LSM_YSZ domain
    model_symbol = eval(Symbol(physical_model_name))
        
    
    if length(physical_model_name) > 13 && physical_model_name[1:13]=="YSZ_LSM_model"             
      push!(prms_names_in, "T")
      push!(prms_values_in, T)
      
      push!(prms_names_in, "pO2")
      push!(prms_values_in, pO2)
      
      push!(prms_names_in, "S_ellyt")
      push!(prms_values_in, AreaEllyt)
      
      #push!(prms_names_in, "h_min")
      #push!(prms_values_in, (10.0^dx_exp)*1.0e-4                                           )
      
      params_dict=Dict{Symbol, Any}(Pair.(Symbol.(prms_names_in), prms_values_in))
      
      sys = model_symbol.sys(params_dict=params_dict)
      #sys = model_symbol.sys()      
      
      #model_symbol.test_IV(sys)
      #model_symbol.impedance_sweep_test(sys)
            
      eq_solution = model_symbol.equilibrium_solution(sys)            
      
      # # # # # # # # # # I-V curve # # # # # # # # # # # # #
# #       begin
# #         bias_range = collect(0.0 : 0.05 : 0.10)
# #         
# #           #bias_range = [0.05]
# #         biased_steadystate_solution = model_symbol.phi_stationary_sweep(sys, eq_solution, bias_range=bias_range)[!, :solution]
# #         for (i, bias) in enumerate(bias_range)
# #           #for i in 1:length(bias_range)  
# #           
# #           model_symbol.plotsolution(sys, biased_steadystate_solution[i], zoom=5.0e-9)
# #           println(" ---------------- ")
# #           println("I-YSZ    = ",model_symbol.YSZ_current_neg(sys, biased_steadystate_solution[i]))
# #           println("I-LSM    = ",model_symbol.LSM_current(sys, biased_steadystate_solution[i]))
# #           println("I-legacy = ",model_symbol.legacy_current(sys, biased_steadystate_solution[i]))
# #           println("I_OXIDE  = ",model_symbol.get_oxide_electric_flux(sys, biased_steadystate_solution[i]))
# #           println("I_elect  = ",model_symbol.get_electron_electric_flux(sys, biased_steadystate_solution[i]))
# #           
# #           pause(1)
# #         end
# #       end

      function evaluate_total_current(generic_current_funcion, U_old, U_new, t_step)
        (C_old_steady, C_old_trans) = generic_current_funcion(sys, U_old)
        (C_new_steady, C_new_trans) = generic_current_funcion(sys, U_new)
          
        return (C_old_steady + C_new_steady)/2 + (C_new_trans - C_old_trans)/t_step
      end
      
      if EIS_IS
        if bias == 0
          biased_steadystate_solution = eq_solution 
          @show biased_steadystate_solution
        else
          df = model_symbol.phi_stationary_sweep(sys, eq_solution, bias_range=collect(0.05 : sign(bias)*0.05 : bias))          
          @show df
          biased_steadystate_solution = last(df)[:solution]          
        end
      
        
        # # # # # # steadystate test # # # # # # #
#         begin      
#           println("YSZ_new  = ",model_symbol.YSZ_testing_current(sys, biased_steadystate_solution))    
#           println("LSM_new  = ",model_symbol.LSM_testing_current(sys, biased_steadystate_solution))              
#                                     
#           println("I-YSZ    = ",model_symbol.YSZ_current_neg(sys, biased_steadystate_solution))
#           println("I-LSM    = ",model_symbol.LSM_current(sys, biased_steadystate_solution))
#           println("I-legacy = ",model_symbol.legacy_current(sys, biased_steadystate_solution))
#           println("I_OXIDE  = ",model_symbol.get_oxide_electric_flux(sys, biased_steadystate_solution))
#           println("I_elect  = ",model_symbol.get_electron_electric_flux(sys, biased_steadystate_solution))          
#         end
        
        
        # # # # # potential step response # # # # # # #
        begin
          tstep = 1e-8
          df_sol = model_symbol.potential_step_response(sys, tstep = tstep)[!, :solution]
          @show length(df_sol)
          
          #for i in collect(3 : 1 : length(df_sol))
          for i in collect(20 : 10 : 60)
            model_symbol.plotsolution(sys, df_sol[i], zoom=5.0e-9)
            
            println("I-YSZ          = ", evaluate_total_current(model_symbol.YSZ_current_neg, df_sol[i-1], df_sol[i], tstep))
            println("I-LSM          = ", evaluate_total_current(model_symbol.LSM_current, df_sol[i-1], df_sol[i], tstep))
            println("I-YSZ_testing  = ", evaluate_total_current(model_symbol.YSZ_testing_current, df_sol[i-1], df_sol[i], tstep))
            println("I-LSM_testing  = ", evaluate_total_current(model_symbol.LSM_testing_current, df_sol[i-1], df_sol[i], tstep))
            #println("I-legacy = ",model_symbol.legacy_current(sys, df_sol[1]))
            pause(5)
          end
          return
        end        

#         # # # # # # impedance test # # # # # # # # # 
#         begin
#           p = Plots.plot(ratio=:equal)
#           for cF in [ 
#                       model_symbol.LSM_current,
#                       model_symbol.YSZ_current_neg
#                       ] 
#                       #model_symbol.YSZ_testing_current, model_symbol.LSM_testing_current]
#               df = model_symbol.impedance_sweep(sys, biased_steadystate_solution,                                 
#                                     currentF=cF,                                                                    
#                                     #currentF=LSM_current,                                              
#                                     )
#               Plots.plot!(p, real.(df.Z), -imag.(df.Z), seriestype=:scatter, label=string(Symbol(cF)))
#           end
#           gui(p)
#           return
#         end
        
      end
      if voltammetry && fast_CV_mode
        
        df = model_symbol.test_IV(sys, bound=upp_bound_eta, step=upp_bound_eta/sample, plot_bool=true,
                                  cF_list=[ model_symbol.LSM_current, 
                                            model_symbol.LSM_testing_current, 
                                            model_symbol.YSZ_testing_current,
                                            model_symbol.YSZ_current_neg]
                    )
        current_col = deepcopy(df[!, :YSZ_testing_current])
        
        
        return DataFrame(U = [], I = []) ##### TODO !!! ! ! !
      end
    end
    
    ##################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if physical_model_name=="ysz_model_GAS_LoMA_generic"      
      extended_LSM_domain_mode=true
    else
      extended_LSM_domain_mode = false
    end        
    
    
    bulk_species = model_symbol.bulk_species
    surface_species = model_symbol.surface_species
    surface_names = model_symbol.surface_names
    iphi = model_symbol.iphi
    iy = model_symbol.iy
    iyAs = model_symbol.iyAs
    iyOs = model_symbol.iyOs
#     if extended_LSM_domain_mode
#       iphiYSZ = model_symbol.iphiYSZ
#     end
    index_driving_species = model_symbol.index_driving_species
    
    
    
    function plot_solution(U, X, x_factor=10^9, x_label="", plotted_length= 5.0; marker="o")
      point_marker_size = 5
      
      
      PyPlot.subplot(211)
      PyPlot.plot((x_factor)*X[:],U[iy,:],label="y")
      for (i, idx) in enumerate(surface_species[1:2])
        PyPlot.plot(0,U[idx,1],marker, markersize=point_marker_size, label=surface_names[i])
      end
      
      PyPlot.ylim(-0.5,1.1)
      PyPlot.xlim(-0.01*plotted_length, plotted_length)
      PyPlot.xlabel(x_label)
      PyPlot.legend(loc="best")
      PyPlot.grid()
      
      PyPlot.subplot(212)
      if extended_LSM_domain_mode
        for i in [3]
          PyPlot.plot(0,U[surface_species[i],1],marker, markersize=point_marker_size, label=surface_names[i])
        end
      end
      PyPlot.plot((x_factor)*X[:],U[iphi,:],label="phi (V)")
      PyPlot.xlim(-0.01*plotted_length, plotted_length)
      PyPlot.xlabel(x_label)
      PyPlot.legend(loc="best")
      PyPlot.grid()
    end
    
    function set_eta(eta)      
      if (
          physical_model_name == "ysz_model_GAS_LoMA" ||
          physical_model_name == "ysz_model_GAS_LoMA_shared" ||
          physical_model_name == "ysz_model_GAS_LoMA_Temperature" ||
          physical_model_name == "ysz_model_GAS_LoMA_generic"
          )
        sys.boundary_values[index_driving_species, 1] = parameters.phi_eq + eta
      else
        model_symbol.set_eta!(sys, eta)
        
        #sys.physics.data.bias = eta
      end
#       if extended_LSM_domain_mode
#         # E^LSM =  zeta*E^YSZ
#         # but more important is ->     index_driving_species = iphiLSM
#         
#         sys.boundary_values[index_driving_species, 1] = parameters.phiLSM_eq + eta*ZETA_fac   
#       else
#         # E^LSM =  zeta* ( E^YSZ + IR )
#         # but more important is ->     index_driving_species = iphiS
#         
#         sys.boundary_values[index_driving_species, 1] = parameters.phi_eq + eta
#       end
    end
    
    function ramp_to_phiS!(U, U0, phiS, phiS0;                          
                          verbose=false)
      got_exception = Nothing
      ramp_isok = false
      
      if verbose
        @show phiS0, phiS
      end
      
     
      
      given_damp_growth = control.damp_growth      
      ramp_nodes_list =       [0,                 0,    0,  4,   4,   16,  32 ]
      damp_growth_in_round =  [given_damp_growth, 1.5, 1.1, 1.5, 1.1, 1.1, 1.1]
            
      phi_ramp_solved = phiS0
      
      U_ramp = deepcopy(U0)
      
      for (round, ramp_nodes) in enumerate(ramp_nodes_list)
        
        control.damp_growth=damp_growth_in_round[round]        
        actual_ramp_step = ( (phiS - phiS0) > 1e-8 ?  (phiS - phiS0) : 1.0 )
        for phi_ramp in (ramp_nodes == 0 ? phiS : collect(phi_ramp_solved : actual_ramp_step/ramp_nodes : phiS))
            try              
              if verbose
                @show phi_ramp                
              end
              set_eta(phi_ramp - parameters.phi_eq)                  
              solve!(U, U_ramp, sys, control=control)
              U_ramp .= U
              phi_ramp_solved = phi_ramp
              
              ramp_isok = true
              
            catch e
              if e isa InterruptException
                rethrow(e)
              else
                #println("fail")
                got_exception = e
                ramp_isok=false
                break
              end
            end
        end
        if ramp_isok
          control.damp_growth = given_damp_growth
          break
        end
      end
      if !(ramp_isok)
        println("ERROR: run_new: ramp cannot reach the steady state") 
        println(got_exception)
        throw(Exception)
      end                    
    end
    
    ############################################
    
    # prms_in = [ A0, R0, DGA, DGR, beta, A ]
    
    #width_Ellyt = 0.00045           # m     width of the half-cell
    #width_Ellyt = 0.0005           # m     width of the half-cell
#     if dlcap
#         AreaEllyt = 1.0      # m^2    
#         if print_bool
#             println("dlcap > area = 1")
#         end
#     end
    #
    dx_start = 10^convert(Float64,dx_exp)
    X=width*VoronoiFVM.geomspace(0.0,1.0,dx_start,1e-1)
    #println("X = ",X)
    #
    grid=VoronoiFVM.Grid(X)
    #
    
    # Physical condtions and parameters setting
    parameters=model_symbol.YSZParameters()
    

    # experimental setting
    parameters.pO2 = pO2
    parameters.T = T
    
    #
    model_symbol.set_parameters!(parameters, prms_values_in, prms_names_in)


    
    # update the "computed" values in parameters
    parameters = model_symbol.YSZParameters_update!(parameters)

    
    #@show parameters.yB
    #@show parameters.DD
  
    # conductivity stuff
    if conductivity_fitting
      sigma = model_symbol.get_conductivity(parameters)
      R_ohm = (1/sigma)*(width/AreaEllyt)
      
      ###
      #
      ######
      #@show sigma
      #@show R_ohm
      
      
      
      return sigma, R_ohm
    end

#     if debug_print_bool 
#         println("NEW ---- ")
#         println("prms = ",prms_in)
#         println("add_prms = (",parameters.DD,",",parameters.nu,",",parameters.nus,",",parameters.ms_par,")")
#         println("rest_prms = (",parameters.T,",",parameters.pO2,")")
#     end
    
    #model_symbol.printfields(parameters)
    
    
    
    # assembling of the system
    if extended_LSM_domain_mode
      physics=VoronoiFVM.Physics(
          data=parameters,
          num_species=size(bulk_species,1)+size(surface_species,1),
          storage=model_symbol.storage!,
          flux=model_symbol.flux!,
          reaction=model_symbol.reaction!,
          breaction=model_symbol.breaction!,
          generic=model_symbol.generic_operator!,
          generic_sparsity=model_symbol.generic_operator_sparsity,
          bstorage=model_symbol.bstorage!
      )
    else
      physics=VoronoiFVM.Physics(
          data=parameters,
          num_species=size(bulk_species,1)+size(surface_species,1),
          storage=model_symbol.storage!,
          flux=model_symbol.flux!,
          reaction=model_symbol.reaction!,
          breaction=model_symbol.breaction!,          
          bstorage=model_symbol.bstorage!
      )    
    end
    
    
    
    
    if print_bool
        model_symbol.printfields(parameters)
    end

    #sys=VoronoiFVM.SparseSystem(grid,physics)
    sys=VoronoiFVM.DenseSystem(grid,physics)
    for idx in bulk_species
      enable_species!(sys, idx, [1])
    end
    for idx in surface_species
      enable_boundary_species!(sys, idx, [1])
    end

    # set boundary conditions
    model_symbol.set_typical_boundary_conditions!(sys, parameters)
    
    
    # get initial value of type unknows(sys) and initial voltage
    inival = model_symbol.get_typical_initial_conditions(sys, parameters)
    phi0 = parameters.phi_eq

    if extended_LSM_domain_mode
      # TODO --- change
      #ZETA_fac = 1
      ZETA_fac = (1 + parameters.e_fac)
    else
      ZETA_fac = (1 + parameters.e_fac)
    end
    

    #
    control=VoronoiFVM.NewtonControl()
    control.verbose=false
    control.tol_linear=1.0e-5
    control.tol_relative=1.0e-8
    #control.tol_absolute=1.0e-4
    control.max_iterations=300
    control.max_lureuse=0
    control.damp_initial=1.0e-3
    control.damp_growth=2.0
    
    
##### used for diplaying initial conditions vs steady state    
#     figure(111)
#     plot_solution(inival, X, 10^9, marker="o")
#     
# #     @show sys.boundary_factors[1]
# #     @show sys.boundary_values[1]
# #     #@show sys.boundary_factors[1] = VoronoiFVM.Dirichlet
# #     #@show sys.boundary_values[1] = parameters.phiS_eq
# #     
# #     
# #     @show inival[iphi,1] 
# #     #@show inival[iphiYSZ,1]
#      
#     
#     
#     steadystate = unknowns(sys)
#     solve!(steadystate, inival, sys, control=control)
#     plot_solution(steadystate, X, 10^9, marker="x")
#     println(sum(steadystate))
#     
#     #@show steadystate[iphi,1]
#     
#     #@show steadystate[iphiYSZ,1]
#     #@show parameters.phiLSM_eq  - steadystate[iphi,1]*(1 + parameters.e_fac)
#     
#     #value_phiLSM = parameters.phiLSM_eq
#     U1 = unknowns(sys)
#     
#     eta_sim = -0.8
#     if extended_LSM_domain_mode
#       sys.physics.data.phiLSM = parameters.phiLSM_eq + eta_sim
#     else
#       sys.boundary_values[index_driving_species, 1]=parameters.phi_eq + eta_sim/ZETA_fac
#     end
#     
#     
#     solve!(U1, steadystate, sys, control=control)
#     plot_solution(U1, X, 10^9, marker="x")
#     
#     if extended_LSM_domain_mode
#       @show steadystate[iphi,1]    
#       @show steadystate[iphiYSZ,1]
#       @show parameters.phiLSM_eq + 1.0  - U1[iphi,1]*(1 + parameters.e_fac)
#     end
#       
#     return
################

    #@show parameters.phi_eq
    #@show parameters.phiS_eq
    #@show parameters.phiLSM_eq
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################################
    ############################################
    ############################################
    ############################################
    ############################################
    ######### dlcap_analytical computation #####

    if dlcap_analytical
      @show parameters.phi_eq
      phi_range = parameters.phi_eq .+ checknodes
      cbl, cs = AreaEllyt.*model_symbol.direct_capacitance(parameters, phi_range)
      out_df = DataFrame(U = phi_range .- parameters.phi_eq, C = cbl + cs, Cs = cs, Cb = cbl) 
      return out_df
    end
    
    #####################################
    #####################################
    #####################################
    #####################################
    #####################################
    ######### EIS computation #############

    if EIS_IS || EIS_TDS
        #
        # Transient part of measurement functional 
        #
        function meas_tran(meas, u)
          U=reshape(u,sys)
          model_symbol.set_meas_and_get_tran_I_contributions(meas, U, sys, parameters, AreaEllyt, X)
        end
        #
        # Steady part of measurement functional
        #
        function meas_stdy(meas, u)
          U=reshape(u,sys)
          model_symbol.set_meas_and_get_stdy_I_contributions(meas, U, sys, parameters, AreaEllyt, X)
        end
        
        #
        # The overall measurement (in the time domain) is meas_stdy(u)+ d/dt meas_tran(u)
        #

        # Calculate steady state solution
        steadystate = unknowns(sys)
        phi_steady = parameters.phi_eq + bias/ZETA_fac      
        
        excited_spec=index_driving_species
        excited_bc=1
        excited_bcval=phi_steady
        
        # relaxation (ramp to phi_steady)
        ramp_to_phiS!(steadystate, inival, phi_steady, parameters.phi_eq)

        
        # Create impedance system
        isys=VoronoiFVM.ImpedanceSystem(sys,steadystate,excited_spec, excited_bc)

        # Derivatives of measurement functionals
        # For the Julia magic behind this we need the measurement functionals
        # as mutating functions writing on vectors.
        dmeas_stdy=measurement_derivative(sys,meas_stdy,steadystate)
        dmeas_tran=measurement_derivative(sys,meas_tran,steadystate)


        
        # Impedance arrays
        z_timedomain=zeros(Complex{Float64},0)
        z_freqdomain=zeros(Complex{Float64},0)
        all_w=zeros(0)

        w = 2*pi*f_range[1]

        
        # Frequency loop
        while w < 2*pi*f_range[2]
            
            print_bool && @show w
            #@show w
            push!(all_w,w)
            if EIS_IS
                # Here, we use the derivatives of the measurement functional
                zfreq=freqdomain_impedance(isys,w,steadystate,excited_spec,excited_bc,excited_bcval,dmeas_stdy, dmeas_tran)
                inductance = im*parameters.L*w
                push!(z_freqdomain, inductance + (1.0/zfreq) )#   /(1 + parameters.e_fac))
                print_bool && @show zfreq
            end
            
            if EIS_TDS
                # Similar API, but use the the measurement functional themselves    
                
                ztime=timedomain_impedance(sys,w,steadystate,excited_spec,excited_bc,excited_bcval,meas_stdy, meas_tran,
                                    tref=tref, excitation_amplitude=1.0e-3,
                                    fit=true, tol_amplitude=1.0e-3, fit_window_size=20.0, plot_amplitude=false)
                inductance = im*parameters.L*w                        
                push!(z_timedomain, inductance + (1.0/ztime))#/(1 + parameters.e_fac))
                print_bool && @show ztime
                @show w, z_timedomain[end]
            end

            
            
            
            # growth factor such that there are 10 points in every order of magnitude
            # (which is consistent with "freq" list below)
            #w=w*1.25892     
            
            w = w*f_range[3]
            
        end
# # # # ##############################
# # # #         pyplot = true 
# # # #         ##########################################
        if pyplot
            EIS_TDS_figure_num = 323
            
            function positive_angle(z)
                ϕ=angle(z)
                if ϕ<0.0
                    ϕ=ϕ+2*π
                end
                return ϕ
            end

            PyPlot.clf()
            PyPlot.subplot(311)
            PyPlot.grid()

            if EIS_IS
                PyPlot.semilogx(all_w,positive_angle.(1.0/z_freqdomain)',label="\$i\\omega\$",color=:red)
            end
            if EIS_TDS
                PyPlot.semilogx(all_w,positive_angle.(1.0/z_timedomain)',label="\$\\frac{d}{dt}\$",color=:green)
            end
            PyPlot.xlabel("\$\\omega\$")
            PyPlot.ylabel("\$\\phi\$")
            PyPlot.legend(loc="upper left")


            PyPlot.subplot(312)
            PyPlot.grid()
            if EIS_IS
                PyPlot.loglog(all_w,abs.(1.0/z_freqdomain)',label="\$i\\omega\$",color=:red)
            end
            if EIS_TDS
                PyPlot.loglog(all_w,abs.(1.0/z_timedomain)',label="\$\\frac{d}{dt}\$",color=:green)
            end
            PyPlot.xlabel("\$\\omega\$")
            PyPlot.ylabel("a")
            PyPlot.legend(loc="lower left")

            
            ax=PyPlot.subplot(313)
            PyPlot.grid()
            ax.set_aspect(aspect=1.0)
            if EIS_IS
                PyPlot.plot(real(z_freqdomain),-imag(z_freqdomain),label="\$i\\omega\$", color=:red)
            end
            if EIS_TDS
                PyPlot.plot(real(z_timedomain),-imag(z_timedomain),label="\$\\frac{d}{dt}\$", color=:green)
            end
            PyPlot.xlabel("Re")
            PyPlot.ylabel("Im")
            PyPlot.legend(loc="lower center")
            PyPlot.tight_layout()
            pause(1.0e-10)                      
        end
        
        if out_df_bool
          if EIS_IS
            EIS_df = DataFrame(f = all_w/(2*pi), Z = z_freqdomain)
            return EIS_df
          elseif EIS_TDS
            EIS_df = DataFrame(f = all_w/(2*pi), Z = z_timedomain)
            return EIS_df
          end
        end
        
# # # # #         ###########
# # # # #         pause(1000)
# # # # #         ##########
    end

    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    
    # code for performing STEP
    
    if STEP
      return "implement zeta e_fac handeling"
# # # #       phi = parameters.phi_eq + bias
# # # #       U = unknowns(sys)
# # # #       U0 = unknowns(sys)
# # # #       only_formal_meas = [0.0,0.1,0.1,0.1,0.1]
# # # #       
# # # #       # solving stationary state
# # # #       sys.boundary_values[index_driving_species,1]=phi
# # # #       solve!(U0,inival,sys, control=control)
# # # # 
# # # #       I_contributions_tran_0 = model_symbol.set_meas_and_get_tran_I_contributions(only_formal_meas, U0, sys, parameters, AreaEllyt, X)
# # # #       I_contributions_stdy_0 = model_symbol.set_meas_and_get_stdy_I_contributions(only_formal_meas, U0, sys, parameters, AreaEllyt, X)
# # # #       #plot_solution(U0, X)
# # # #       #pause(2)
# # # #       
# # # #       # relaxation
# # # #       total_charge = 0
# # # #       
# # # #       t = 0
# # # #       dt = dt_start
# # # #       phi = parameters.phi_eq
# # # #       sys.boundary_values[index_driving_species,1]=phi
# # # #       while t < t_end
# # # #         t += dt
# # # #         control.Δt = dt
# # # #         evolve!(U, U0, sys, [0, dt], control=control)
# # # # 
# # # #         #@show t
# # # #         #plot_solution(U, X)
# # # #         
# # # #         I_contributions_tran = model_symbol.set_meas_and_get_tran_I_contributions(only_formal_meas, U, sys, parameters, AreaEllyt, X)
# # # #         I_contributions_stdy = model_symbol.set_meas_and_get_stdy_I_contributions(only_formal_meas, U, sys, parameters, AreaEllyt, X)
# # # #         I = (
# # # #           sum((I_contributions_stdy .+ I_contributions_stdy_0)./2)
# # # #           +
# # # #           sum((I_contributions_tran .- I_contributions_tran_0)./dt)
# # # #         )
# # # #         
# # # #         total_charge += I*dt
# # # #         @show I*dt
# # # #         
# # # #         I_contributions_tran_0 = I_contributions_tran
# # # #         I_contributions_stdy_0 = I_contributions_stdy
# # # #         U0 .= U
# # # #         dt *= dt_fac
# # # #       end
# # # #       return total_charge
    end
    
    
    
    
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    
    # code for performing CV(f) ---- FAST cv MODE
    
    if voltammetry && fast_CV_mode   
        upp_bound_phi = upp_bound_eta/ZETA_fac
        low_bound_phi = low_bound_eta/ZETA_fac
        
        tstep=Float64(((upp_bound_phi-low_bound_phi)/2)/voltrate/sample)
        
        if save_files || out_df_bool          
            out_df = DataFrame(t = Float64[], U = Float64[], I = Float64[], Ir = Float64[], Is = Float64[], Ib = Float64[], Ibb = Float64[])          
        end
  
        only_formal_meas = [0.0]  # to pass as an argument to set_meas_get_tran_I_contributions

        U = unknowns(sys)        
        U0 = deepcopy(inival)              
              
        # calculating directly steadystate ###############                                        
        set_eta(0.0)
        solve!(U0,inival,sys, control=control)
        U_eq = deepcopy(U0)

        
        if save_files || out_df_bool
          push!(out_df, (0, 0, 0, 0, 0, 0, 0))
        end            
        
        # compute the positive half
        istep = 1
        dir = 1
        phi_prev = parameters.phi_eq
        phi = parameters.phi_eq + voltrate*tstep*dir        
        while phi <= (upp_bound_phi + parameters.phi_eq) + upp_bound_phi/100.
          
          ramp_to_phiS!(U, U0, phi, phi_prev)
          
          # Steady part of measurement functional        
          I_contributions_stdy = model_symbol.set_meas_and_get_stdy_I_contributions(only_formal_meas, U, sys, parameters, AreaEllyt, X)
          #
          Ir= I_contributions_stdy[1]          
          if save_files || out_df_bool
            push!(out_df,[istep*tstep,   ZETA_fac*(phi - parameters.phi_eq),   Ir+0,  Ir,  0,  0,  0])
          end
          
          # preparing for the next step
          istep+=1
          U0.=U
          phi_prev = phi
          phi+=voltrate*tstep*dir       
        end
        
        # complete the other direction
        for i in collect(size(out_df)[1]: -1 : 1)        
          push!(out_df, (istep*tstep, Tuple(out_df[i, 2:end])...))          
          istep+=1
        end        
        
        
        # compute the negative half
        istep_negative = istep
        U0 .= U_eq
        dir = -1        
        phi_prev = parameters.phi_eq
        phi = parameters.phi_eq + voltrate*tstep*dir
        while phi >= (low_bound_phi + parameters.phi_eq) - upp_bound_phi/100.
          
          ramp_to_phiS!(U, U0, phi, phi_prev)
          
          # Steady part of measurement functional        
          I_contributions_stdy = model_symbol.set_meas_and_get_stdy_I_contributions(only_formal_meas, U, sys, parameters, AreaEllyt, X)
          #
          Ir= I_contributions_stdy[1]          
          if save_files || out_df_bool
            push!(out_df,[istep*tstep,   ZETA_fac*(phi - parameters.phi_eq),   Ir+0,  Ir,  0,  0,  0])
          end
          
          # preparing for the next step
          istep+=1
          U0.=U
          phi_prev = phi
          phi+=voltrate*tstep*dir     
        end

        # complete the other direction 
        for i in collect(size(out_df)[1] : -1 : istep_negative)   
          push!(out_df, (istep*tstep, Tuple(out_df[i, 2:end])...))          
          istep+=1
        end 
        
        if out_df_bool
          return out_df
        else
          return
        end      
    end
    
    
    
    
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################
    
    
    
    # code for performing CV    
    if voltammetry
        upp_bound_phi = upp_bound_eta/ZETA_fac
        low_bound_phi = low_bound_eta/ZETA_fac
    
        istep=1
        phi=0
        phi_prev = 0
        phi_out = 0

        # inicializing storage lists
        y0_range=zeros(0)
        surface_species_range=Array{Float64}(undef, (size(surface_species,1), 0))
        yOs_range=zeros(0)
        phi_range=zeros(0)
        #
        Is_range=zeros(0)
        Ib_range=zeros(0)
        Ibb_range=zeros(0)
        Ir_range=zeros(0)
        
        
        if save_files || out_df_bool
          if dlcap
            out_df = DataFrame(t = Float64[], U = Float64[], C = Float64[], Cr = Float64[], Cs = Float64[], Cb = Float64[], Cbb = Float64[])
          else
            out_df = DataFrame(t = Float64[], U = Float64[], I = Float64[], Ir = Float64[], Is = Float64[], Ib = Float64[], Ibb = Float64[])
          end
        end
        
        cv_cycles = 1
        sinus_mode = false
        
        relax_counter = 0
        istep_cv_start = -1
        time_range = zeros(0)  # [s]
        
        only_formal_meas = [0.0]  # to pass as an argument to set_meas_get_tran_I_contributions

        if print_bool
            print("calculating linear potential sweep\n")
        end
        direction_switch_count = 0

        tstep=Float64(((upp_bound_phi-low_bound_phi)/2)/voltrate/sample)
        if print_bool
            @printf("tstep %g = \n", tstep)
        end
        if phi0 > 0
            dir=1
        else
            dir=-1
        end
        
        if pyplot
            PyPlot.close()
            PyPlot.ion()            
            PyPlot.figure(321, figsize=(10,8))
            #PyPlot.figure(figsize=(5,5))
        end
        
        
        U = unknowns(sys)
        U0 = unknowns(sys)
        U0 .=inival
        if test
            U .= inival
        end
        
              
        if true
          # calculating directly steadystate ###############                  
          set_eta(0.0 - parameters.phi_eq)
          solve!(U0,inival,sys, control=control)
          if save_files || out_df_bool
            push!(out_df, (0, 0, 0, 0, 0, 0, 0))
          end
        
          append!(y0_range,U0[iy,1])
          surface_species_to_append = []
          for (i, idx) in enumerate(surface_species)
            append!(surface_species_to_append,U0[idx,1])
          end
          surface_species_range = hcat(surface_species_range, surface_species_to_append)
          
          append!(phi_range, U0[iphi, 1])
          #
          I_contributions_stdy_0 = [0]
          I_contributions_tran_0 = [0, 0, 0]
          append!(Ib_range,0)
          append!(Is_range,0)
          append!(Ibb_range,0)
          append!(Ir_range,0)
          #
          append!(time_range,0)
          
          
          # preparing fist step          
          istep_cv_start = 0          
          dir=1
          state = "relaxation"
          phi=parameters.phi_eq
          relaxation_length = 2*2.0/sample
          
          phi_prev = parameters.phi_eq
          #phi = phi + voltrate*dir*tstep
          phi_out = (phi_prev + phi)/2.0
          ##################################################
          # using the ramp ... not recommended
        else
          state = "ramp"
          I_contributions_stdy_0 = [0]
          I_contributions_tran_0 = [0, 0, 0]
          if print_bool
              println("phi_equilibrium = ",phi0)
              println("ramp ......")
          end
        end
        

        
        while state != "cv_is_off"                           
            if state=="ramp" && ((dir==1 && phi > phi0) || (dir==-1 && phi < phi0))
                phi = phi0
                state = "relaxation"
                if print_bool
                    println("relaxation ... ")
                end
            end  
            if state=="relaxation" && relax_counter >= sample*relaxation_length
                relax_counter += 1
                if sinus_mode
                  state="sinus_mode"
                else
                  state="cv_is_on"
                end
                istep_cv_start = istep
                dir=1
                if print_bool
                    print("cv ~~~ direction switch: ")
                end
            end
            if state=="sinus_mode" && istep > 100
              state = "cv_is_off"
            end
            if state=="cv_is_on" && (phi <= low_bound_phi-0.00000001+phi0 || phi >= upp_bound_phi+0.00000001+phi0)
                dir*=(-1)
                # uncomment next line if phi should NOT go slightly beyond limits
                #phi+=2*voltrate*dir*tstep
            
                direction_switch_count +=1
                if print_bool
                    print(direction_switch_count,", ")
                end
            end       
            if state=="cv_is_on" && (dir > 0) && (phi > phi0 + voltrate*dir*tstep + 0.000001) && (direction_switch_count >=2*cv_cycles)
                state = "cv_is_off"
            end
            
            # tstep to potential phi
            set_eta(phi - parameters.phi_eq)
                        
            control.Δt = tstep

            evolve!(U, U0, sys, [0, tstep], control=control)
            
#     #############
#     #TODO EREASE
#     @show U[iphiYSZ,1]
#     @show U[iphi,1]
#     @show U[iphiLSM,1]    
#     #############
        
            
            # Transient part of measurement functional
            I_contributions_tran = model_symbol.set_meas_and_get_tran_I_contributions(only_formal_meas, U, sys, parameters, AreaEllyt, X)
            #
            # Steady part of measurement functional
            #
            I_contributions_stdy = model_symbol.set_meas_and_get_stdy_I_contributions(only_formal_meas, U, sys, parameters, AreaEllyt, X)
            
            # reaction average
            Ir= 0.5*(I_contributions_stdy[1] + I_contributions_stdy_0[1])
            
            # time derivatives           
            Is  = (I_contributions_tran[1] - I_contributions_tran_0[1])/tstep                
            Ib  = (I_contributions_tran[2] - I_contributions_tran_0[2])/tstep 
            Ibb = (I_contributions_tran[3] - I_contributions_tran_0[3])/tstep
            
          
          
            #@show state
            #@show Ir+Is+Ib+Ibb
            if false && !(test || test_from_above)
                @printf("t = %.2g   U = %g   ys0 = %g  yAs = %g  yOs = %g  r=(%g, %g, %g)  I=(%g, %g, %g, %g)\n", istep*tstep, phi, U[iy,1], U[iyAs,1], U[iyOs,1], 
                      model_symbol.exponential_oxide_adsorption(parameters, U[:,1], debug_bool=true),
                      model_symbol.electroreaction(parameters, U[:,1], debug_bool=true),
                      model_symbol.exponential_gas_adsorption(parameters, U[:,1], debug_bool=true),
                      
                      Ir, Is, Ib, Ibb
                      )
            end
            
            # storing data
            append!(y0_range,U[iy,1])
            
            surface_species_to_append = []
            for (i, idx) in enumerate(surface_species)
              append!(surface_species_to_append,U[idx,1])
            end
            surface_species_range = hcat(surface_species_range, surface_species_to_append)
            
            append!(phi_range,U[iphi,1])              
            append!(time_range,tstep*istep)
            if istep > 2              
              append!(Ib_range,Ib)
              append!(Is_range,Is)
              append!(Ibb_range,Ibb)
              append!(Ir_range, Ir)              
            else              
              append!(Ib_range,0)
              append!(Is_range,0)
              append!(Ibb_range,0)
              append!(Ir_range, 0)                            
            end
            
            if state=="cv_is_on"
                if save_files || out_df_bool
                    if dlcap
                        push!(out_df,[(istep-istep_cv_start)*tstep, phi_out-phi0, 
                          dir*(Ir+Is+Ib+Ibb)/voltrate, 
                          dir*Ir/voltrate, 
                          dir*Is/voltrate,
                          dir*Ib/voltrate,                            
                          dir*Ibb/voltrate])
                    else
                        push!(out_df,[(istep-istep_cv_start)*tstep,   ZETA_fac*(phi_out-phi0),   Ib+Is+Ir+Ibb,  Ir,  Is,  Ib,  Ibb])
                    end
                end
            end
                      
            # plotting                  
            if pyplot && istep%2 == 0
                
                figure(321)
                num_subplots=4
                ys_marker_size=6
                PyPlot.subplots_adjust(hspace=0.5)
            
                PyPlot.clf() 
                
                if num_subplots > 0

                      subplot(num_subplots*100 + 11)
                      plot((10^9)*X[:],U[iphi,:],label="phi (V)")
                      plot((10^9)*X[:],U[iy,:],label="y")
                      for (i, idx) in enumerate(surface_species)
                        plot(0,U[idx,1],"o", markersize=ys_marker_size, label=surface_names[i])
                      end
                      l_plot = 5.0
                      PyPlot.xlim(-0.01*l_plot, l_plot)
                      PyPlot.ylim(-0.5,1.1)
                      PyPlot.xlabel("x (nm)")
                      PyPlot.legend(loc="best")
                      PyPlot.grid()
                end
                
                if (long_domain_plot_bool=true)
                    if num_subplots > 1
                        subplot(num_subplots*100 + 12)
                        plot((10^3)*X[:],U[iphi,:],label="phi (V)")
                        plot((10^3)*X[:],U[iy,:],label="y")
                        for (i, idx) in enumerate(surface_species)
                          plot(0,U[idx,1],"o", markersize=ys_marker_size, label=surface_names[i])
                        end
                        PyPlot.ylim(-0.5,1.1)
                        PyPlot.xlabel("x (mm)")
                        PyPlot.legend(loc="best")
                        PyPlot.grid()
                    end
                else
                    if (num_subplots > 1) && (istep_cv_start > -1)
                        cv_range = (istep_cv_start+1):length(phi_range)
                        subplot(num_subplots*100 + 12)
                        plot(phi_range[cv_range].-phi0, ((Is_range + Ib_range + Ir_range + Ibb_range)[cv_range]) ,label="total current")
                        
                        PyPlot.xlabel(L"\eta \ (V)")
                        PyPlot.ylabel(L"I \ (A)")
                        PyPlot.legend(loc="best")
                        PyPlot.grid()
                    end
                end
                
                if num_subplots > 2
                    subplot(num_subplots*100 + 13)
                    plot(time_range,phi_range,label="phi_S (V)")
                    plot(time_range,y0_range,label="y(0)")
                    for (i, idx) in enumerate(surface_species)
                      plot(time_range,surface_species_range[i,:], label=surface_names[i])
                    end
                    PyPlot.xlabel("t (s)")
                    PyPlot.legend(loc="best")
                    PyPlot.grid()
                end
                
                if num_subplots > 3
                    subplot(num_subplots*100 + 14)
                    plot(time_range,Is_range + Ib_range + Ir_range + Ibb_range,label="total I (A)")
                    PyPlot.xlabel("t (s)")
                    PyPlot.legend(loc="best")
                    PyPlot.grid()
                end
                                
                pause(1.0e-10)
            end
            
            # preparing for the next step
            istep+=1
            U0.=U
            I_contributions_stdy_0 = I_contributions_stdy
            I_contributions_tran_0 = I_contributions_tran
            if state=="relaxation"
                relax_counter += 1
                #println("relaxation ... ",relax_counter/sample*100,"%")
            else
                phi_prev = phi
                if sinus_mode
                  phi = parameters.phi_eq + sin(voltrate*tstep*istep)
                else
                  phi+=voltrate*dir*tstep
                end
                phi_out = (phi_prev + phi)/2.0
            end
        end
        
        
        
        # the finall plot
        if pyplot || pyplot_finall
            if pyplot_finall
                PyPlot.figure(figsize=(10,8))
            else
                PyPlot.pause(5)
            end
            
            PyPlot.clf()
            #PyPlot.close()
            #PyPlot.figure(figsize=(5,5))
            
            cv_range = (istep_cv_start+1):length(phi_range)


            subplot(221)
            if dlcap
                plot(phi_range[cv_range].-phi0,( Ib_range[cv_range] )/voltrate,"blue", label="bulk")
                #plot(phi_range[cv_range].-phi0,( Ibb_range[cv_range])/voltrate,label="bulk_grad")
                plot(phi_range[cv_range].-phi0,( Is_range[cv_range] )/voltrate,"green", label="surf")
                plot(phi_range[cv_range].-phi0,( Ir_range[cv_range]  )/voltrate,"red", label="reac")
            else
                plot(phi_range[cv_range].-phi0, Ib_range[cv_range] ,"blue", label="bulk")
                #plot(phi_range[cv_range].-phi0, Ibb_range[cv_range] ,label="bulk_grad")
                plot(phi_range[cv_range].-phi0, Is_range[cv_range] ,"green",label="surf")
                plot(phi_range[cv_range].-phi0, Ir_range[cv_range] ,"red",label="reac")
            end
            if dlcap
                PyPlot.xlabel("nu (V)")
                PyPlot.ylabel(L"Capacitance (F/m$^2$)")  
                PyPlot.legend(loc="best")
                #PyPlot.xlim(-0.5, 0.5)
                #PyPlot.ylim(0, 5)
                PyPlot.grid()
                PyPlot.show()
                PyPlot.pause(5)
                
                PyPlot.clf()
                plot(phi_range[cv_range].-phi0,( (Ib_range+Is_range+Ir_range)[cv_range]  )/voltrate,"brown", label="total")
                PyPlot.xlabel("nu (V)")
                PyPlot.ylabel(L"Capacitance (F/m$^2$)") 
                PyPlot.legend(loc="best")
                #PyPlot.xlim(-0.5, 0.5)
                #PyPlot.ylim(0, 5)
                PyPlot.grid()
                PyPlot.show()
                #PyPlot.pause(10)
            else
                PyPlot.xlabel(L"\eta \ (V)")
                PyPlot.ylabel(L"I \ (A)")
                PyPlot.legend(loc="best")
                PyPlot.grid()
            end
            
            
            subplot(222)
            if dlcap  
                cbl, cs = model_symbol.direct_capacitance(parameters, collect(float(low_bound_phi):0.001:float(upp_bound_phi)))
                plot(collect(float(low_bound_phi):0.001:float(upp_bound_phi)), (cbl+cs), label="tot CG") 
                plot(collect(float(low_bound_phi):0.001:float(upp_bound_phi)), (cbl), label="b CG") 
                plot(collect(float(low_bound_phi):0.001:float(upp_bound_phi)), (cs), label="s CG") 
                # rescaled by voltrate
                plot(phi_range[cv_range].-phi0, ((Is_range + Ib_range + Ir_range + Ibb_range)[cv_range])/voltrate ,label="rescaled total current")
            else
                plot(phi_range[cv_range].-phi0, ((Is_range + Ib_range + Ir_range + Ibb_range)[cv_range]) ,label="total current")
            end
            PyPlot.xlabel(L"\eta \ (V)")
            PyPlot.ylabel(L"I \ (A)")
            PyPlot.legend(loc="best")
            PyPlot.grid()
            
            
            subplot(223)
            #plot(phi_range, Ir_range ,label="spec1")
            plot(time_range,phi_range,label="phi_s (V)")        
            plot(time_range,y0_range,label="y(0)")
            for (i, idx) in enumerate(surface_species)
              plot(time_range,surface_species_range[i,:], label=surface_names[i])
            end
            PyPlot.xlabel("t (s)")
            PyPlot.legend(loc="best")
            PyPlot.grid()
            
            
            subplot(224, facecolor="w")
            height=0.0
            shift=0.0
            swtch=false
            for name in fieldnames(typeof(parameters))
                if (string(name) == "A0" || swtch)
                    swtch = true
                    value = @sprintf("%.6g", parse(Float64,string(getfield(parameters,name))))
                    linestring = @sprintf("%s: %s", name, value)
                    PyPlot.text(0.01+shift, 0.95+height, linestring, fontproperties="monospace")
                    height+=-0.05
                    if string(name) == "e0" 
                        shift+=0.5
                        height=0.0
                    end
                    if string(name) == "A"
                        PyPlot.text(0.01+shift, 0.95+height, " ", fontproperties="monospace")
                        height+=-0.05
                    end
                end
            end
            parn = ["verbose" ,"pyplot", "width", "voltammetry", "voltrate", "low_bound_phi", "upp_bound_phi", "sample", "phi0"]
            parv =[verbose ,pyplot, width, voltammetry, voltrate, low_bound_phi, upp_bound_phi, sample, @sprintf("%.6g",phi0)]
            for ii in 1:length(parn)
                    linestring=string(parn[ii],": ",parv[ii])
                    PyPlot.text(0.01+shift, 0.95+height, linestring, fontproperties="monospace")
                    height+=-0.05
            end
        end


        
        if save_files
            out_name=string(
            "A0",@sprintf("%.0f",prms_in[1]),
            "_R0",@sprintf("%.0f",prms_in[2]),
            "_GA",@sprintf("%.0f",prms_in[3]),
            "_GR",@sprintf("%.0f",prms_in[4]),
            "_be",@sprintf("%.0f",prms_in[5]),
            "_A",@sprintf("%.0f",prms_in[6]),
            "_vrate",@sprintf("%.2g",voltrate),
            )

            out_data_dir = "./results/CV_data/"
            
            if !ispath(out_data_dir)
                mkpath(out_data_dir)
            end

            CSV.write(string(out_data_dir, out_name,".csv"),out_df)
            
            
            if pyplot
                out_fig_dir = "./results/CV_images/"
            
                if !ispath(out_fig_dir)
                    mkpath(out_fig_dir)
                end
            
                PyPlot.savefig(string(out_fig_dir, out_name,".png"))
            end
        end
        if test
            I1 = integrate(sys, model_symbol.reaction!, U)
            #println(I1)
            return I1[1]
        end
        if out_df_bool
            return out_df
        end
    end
end


end
